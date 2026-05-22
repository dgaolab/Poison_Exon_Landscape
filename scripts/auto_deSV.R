suppressPackageStartupMessages({
  library(matrixStats)
  library(future.apply)
  library(MASS)
  library(DESeq2)
  library(sva)
})


is.fullrank = function(modMat){
  qrx = qr(modMat)
  if (qrx$rank == ncol(modMat)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


fit_single_model = function(x, meta_data = NULL, approach = NULL, fit_formula = NULL, offset_str = NULL, regression_caller = NULL, distro = NULL, defaultlink = NULL, control_args = NULL){
  meta_data[["Y"]] = x
  fit_formula_final = paste0("Y", fit_formula)
  res = tryCatch(
    {
      if (is.null(offset_str)){
        if (approach == 0){
          #non-NB distributions without offset called by glm, or non-gaussian and non-NB distribution without offset called by glmer
          regression_caller(formula(fit_formula_final), data = meta_data, family = distro)
        } else if ((approach == 1) || (approach == 2)){
          #NB distribution without offset called by glm.nb, or gaussian distribution without offset called by lmer, or NB distribution without offset called by glmer.nb
          regression_caller(formula(fit_formula_final), data = meta_data)
        }
      } else {
        if (approach == 0){
          #non-NB distributions with offset called by glm, or non-gaussian distribution with offset called by glmer
          regression_caller(formula(fit_formula_final), data = meta_data, family = distro, offset = defaultlink(meta_data[[offset_str]]))
        } else if (approach == 1){
          #Gaussian distribution with offset called by lmer, or NB distribution with offset called by glmer.nb
          regression_caller(formula(fit_formula_final), data = meta_data, offset = defaultlink(meta_data[[offset_str]]))
        } else if (approach == 2){
          #NB distribution without offset called by glm.nb
          regression_caller(formula(paste0(fit_formula_final, " + offset(log(", offset_str, "))")), data = meta_data)
        }
      }
    }, 
    error = function(e){return(NULL)}
  )
  return(res)
}


model_coef_autocomplete = function(single_model, expected_coef_names, extract_column = NULL){
  coef_mat0 = coef(summary(single_model))
  real_coef_names = rownames(coef_mat0)
  for (coef_name in expected_coef_names){
    if (!(coef_name %in% real_coef_names)){
      coef_mat0 = rbind(coef_mat0, c(NA, NA, NA, NA))
      rownames(coef_mat0)[nrow(coef_mat0)] = coef_name
    }
  }
  if (is.null(extract_column)){
    coef_mat1 = coef_mat0[expected_coef_names, ]
  } else if ((extract_column >= 1) && (extract_column <= 4)){
    coef_mat1 = coef_mat0[expected_coef_names, extract_column]
  }
  
  return(coef_mat1)
}


model_vcov_autocomplete = function(single_model, expected_vcov_names){
  vcov_mat0 = vcov(single_model)
  real_vcov_names = rownames(vcov_mat0)
  for (vcov_name in expected_vcov_names){
    if (!(vcov_name %in% real_vcov_names)){
      vcov_mat0 = rbind(vcov_mat0, rep(NA, ncol(vcov_mat0)))
      rownames(vcov_mat0)[nrow(vcov_mat0)] = vcov_name
      vcov_mat0 = cbind(vcov_mat0, rep(NA, nrow(vcov_mat0)))
      colnames(vcov_mat0)[ncol(vcov_mat0)] = vcov_name
    }
  }
  vcov_mat1 = vcov_mat0[expected_vcov_names, expected_vcov_names]
  
  return(vcov_mat1)
}


merge_vcov_from_model_list = function(model_list, expected_vcov_names){
  tt = lapply(model_list, model_vcov_autocomplete, expected_vcov_names = expected_vcov_names)
  tt1 = lapply(names(tt), FUN = function(x){tmp = as.data.frame(as.matrix(tt[[x]])); rownames(tmp) = NULL; tmp[["name"]] = x; return(tmp)}) 
  tt2 = do.call(rbind, tt1)
  return(tt2)
}


#Define a new S4 data structure to store key information from GLM/DESeq2
setClass("GLMres", representation(y = "matrix", colData = "data.frame", design = "formula", beta = "matrix", stderr = "matrix", pval = "matrix", fdr = "matrix", vcov = "data.frame", resid = "matrix", deSV = "matrix", status = "matrix", args = "list", models = "list"))
setMethod(f = "show", signature = "GLMres", definition = function(object){print("class: GLMres");print("slots: y, colData, design, beta, stderr, pval, fdr, vcov, resid, deSV, status, args, models");print(paste0("deSV dim: ", nrow(object@deSV), " ", ncol(object@deSV)))})


buildDEmodel = function(feature_matrix, sampleTable, method = "DESeq2", estimate_offset = TRUE, offset_str = NULL, apply_sva = TRUE, n_sv = NULL, test_terms = NULL, remove_terms = NULL, random_terms = NULL, distro_family = "poisson", save_models = TRUE){
  #te_name = x[1]
  #print(te_name)
  if (all(test_terms %in% colnames(sampleTable)) == FALSE){stop("'test_terms' must be found in the column names of 'sampleTable'.")}
  if (all(remove_terms %in% colnames(sampleTable)) == FALSE){stop("'remove_terms' must be found in the column names of 'sampleTable'.")}
  #if (all(random_terms %in% colnames(sampleTable)) == FALSE){stop("'random_terms' must be found in the column names of 'sampleTable'.")}
  if (all(remove_terms %in% test_terms) == FALSE){stop("'remove_terms' must be a subset of 'test_terms'.")}
  
  irm = as.matrix(feature_matrix)
  irc = as.data.frame(sampleTable)
  fml = paste0("~",paste(c(test_terms), collapse = " + "))
  mod_tmp = model.matrix(formula(fml), data = irc)
  if (is.fullrank(mod_tmp) == FALSE){stop("The model matrix derived from 'test_terms' is not full rank.")}
  sv_terms = c()
  remove_cols = c()

  message(method, " approach selected.")
  if ((method == "deseq2") || (method == "DESeq2")){
    #DESeq2 based approach for gene expression and exon splicing
    message("Constructing DESeq2 object...")
    mode(irm) = "integer"
    ddssva = DESeqDataSetFromMatrix(irm, irc, ~1)
    if (estimate_offset){
      message("Estimating size factors...")
      ddssva = estimateSizeFactors(ddssva)
      offset_str = "buildDEmodel.sizeFactor"
      colData(ddssva)[[offset_str]] = sizeFactors(ddssva)
    } else {
      if (is.null(offset_str)){
        message("No size factor information provided. Setting to 1...")
        offset_str = "buildDEmodel.sizeFactor"
        colData(ddssva)[[offset_str]] = rep(1, nrow(irc))
        sizeFactors(ddssva) = colData(ddssva)[[offset_str]]
      } else {
        if (all(offset_str %in% colnames(sampleTable)) == FALSE){stop("'offset_str' must be found in the column names of 'sampleTable'.")}
        message("Using '", offset_str, "' column in the design matrix as size factors...")
        sizeFactors(ddssva) = colData(ddssva)[[offset_str]]
      }
    }
    
    if (apply_sva){
      message("Identifying surrogate variables...")
      if (is.null(remove_terms)){
        null_model_str = "1"
      } else {
        null_model_str = paste(c(remove_terms), collapse = " + ")
      }
      mod  = model.matrix(formula(fml), colData(ddssva))
      mod0 = model.matrix(formula(paste0("~",null_model_str)), colData(ddssva))
      message("svaseq FULL model: ", formula(fml))
      message("svaseq null model: ", formula(paste0("~",null_model_str)))
      svseq  = svaseq(counts(ddssva, normalized = T), mod, mod0, n_sv)

      for (i in 1:dim(svseq$sv)[2]){ddssva[[paste0("SV",i)]] = svseq$sv[,i]; fml = paste0(fml," + SV",i); sv_terms = c(sv_terms, paste0("SV",i))}
    }
    
    message("Initiating modeling...")
    design(ddssva) = formula(fml)

    message("Fitting ", method, " models with NB distribution: ", fml, " ...")
    ddssva = DESeq(ddssva, modelMatrixType = "standard", parallel = TRUE)
    
    if (length(c(remove_terms, sv_terms)) > 0){
      message("Correcting against unwanted effects: ", paste(c(remove_terms, sv_terms), collapse = ", ")," ...")
      m = model.matrix(formula(fml), data = colData(ddssva))
      m = as.data.frame(m)
      for (keyword in remove_terms){
        remove_cols = c(remove_cols, which(grepl(paste0("^", keyword), colnames(m))))
      }
      for (svword in sv_terms){
        remove_cols = c(remove_cols, which(colnames(m) == svword))
      }
      m1 = m[,remove_cols]
      X = as.matrix(m1)
      beta = coef(ddssva)[, remove_cols]
      beta = as.matrix(beta)
      y = log2(counts(ddssva, normalized = T))
      cleany = y - t(as.matrix(X) %*% t(beta))
    } else {
      y = log2(counts(ddssva, normalized = T))
      cleany = y
    }
    #cleanY = 2^cleany
    
    message("Collecting results...")
    ddssva@metadata$deSV = cleany

    message("Done.")
    return(ddssva)
    
  } else {
    #matrix-based approach
    #message("Removing lowly variable entries...")
    #row_vars = rowVars(irm)
    #irm = irm[which(row_vars > 1e-4),]

    if (estimate_offset){
      message("Estimating size factors...")
      offset_str = "buildDEmodel.sizeFactor"
      irc[[offset_str]] = estimateSizeFactorsForMatrix(irm)
    } else {
      if (is.null(offset_str)){
        message("No size factor information provided. Will NOT use 'offset' option during modeling.")
      } else {
        if (all(offset_str %in% colnames(sampleTable)) == FALSE){stop("'offset_str' must be found in the column names of 'sampleTable'.")}
        message("Using '", offset_str, "' column in the design matrix as size factors...")
      }
    }
    
    if ((distro_family == "poisson") || (distro_family == "Poisson")){
      #if (!all(sapply(irm, `%%`, 1) == 0)){stop(distro_family, " distribution is selected. The 'feature_matrix' must be all integers.")}
      distro_family = "poisson"
      defaultlink = base::log
      reverselink = base::exp
      distro = poisson
      if (is.null(offset_str)){stop(distro_family, " distribution is selected, but no size factor information provided. Either use 'estimate_offset = TRUE' or provide 'offset_str'.")}
    } else if ((distro_family == "nb") || (distro_family == "NB") || (distro_family == "negative.binomial")){
      #if (!all(sapply(irm, `%%`, 1) == 0)){stop(distro_family, " distribution is selected. The 'feature_matrix' must be all integers.")}
      distro_family = "NB"
      defaultlink = base::log
      reverselink = base::exp
      distro = NULL
      if (is.null(offset_str)){stop(distro_family, " distribution is selected, but no size factor information provided. Either use 'estimate_offset = TRUE' or provide 'offset_str'.")}
    } else if ((distro_family == "gaussian") || (distro_family == "Gaussian")){
      #if (all(sapply(irm, `%%`, 1) == 0)){warning(distro_family, " distribution is selected, but the values in 'feature_matrix' are all integers. Consider using Poisson distribution instead?")}
      distro_family = "gaussian"
      defaultlink = identity
      reverselink = identity
      distro = gaussian
      if (!is.null(offset_str)){
        warning(distro_family, " distribution is selected. The provided 'offset_str' option will be ignored.")
        offset_str = NULL
      }
    }

    sv_model_str = ""
    if (apply_sva){
      message("Identifying surrogate variables...")
      if (is.null(remove_terms)){
        null_model_str = "1"
      } else {
        null_model_str = paste(c(remove_terms), collapse = " + ")
      }
      mod  = model.matrix(formula(fml), irc)
      mod0 = model.matrix(formula(paste0("~",null_model_str)), irc)
      message("Getting rescaled data...")
      if ((distro_family == "gaussian") || (distro_family == "Gaussian")){
        irm_norm = irm
      } else {
        irm_norm = t(apply(irm, 1, FUN = function(x){return(x / as.vector(irc[[offset_str]]))}))
      }
      message("svaseq FULL model: ", formula(fml))
      message("svaseq null model: ", formula(paste0("~",null_model_str)))
      #svseq  = svaseq(as.matrix(irm_norm), mod, mod0, n_sv)
      if (is.null(n_sv)){n_sv = num.sv(dat = defaultlink(irm_norm), mod = mod, method = "be", vfilter = NULL)}
      svseq = irwsva.build(dat = defaultlink(irm_norm), mod = mod, mod0 = mod0, n.sv = n_sv, B = 5)

      for (i in 1:dim(svseq$sv)[2]){irc[[paste0("SV",i)]] = svseq$sv[,i]; sv_model_str = paste0(sv_model_str," + SV",i); sv_terms = c(sv_terms, paste0("SV",i))}
    }
    
    message("Initiating modeling...")
    sp = 0
    if ((method == "glm") || (method == "GLM")){
      fml_fit = paste0(fml, sv_model_str)
      fml_fixed = fml_fit

      if ((distro_family == "poisson") || (distro_family == "gaussian")){
        regression_caller = glm
        #control_args = list()
      } else if (distro_family == "NB"){
        regression_caller = glm.nb
        #control_args = glm.control()
        sp = 2
      }
    } else if ((method == "glmm") || (method == "GLMM") || (method == "lme4") || (method == "glmer") || (method == "lmer") || (method == "mixed")){
      if (is.null(random_terms) || (nchar(random_terms) == 0)){stop("'glmm' is selected so that 'random_terms' cannot be empty or NULL.")}
      #if (all(random_terms %in% test_terms) == FALSE){stop("'random_terms' must be a subset of 'test_terms'.")}
      
      main_terms = test_terms[which(!(test_terms %in% random_terms))]
      main_model_str = paste0("~", paste(c(main_terms), collapse = " + "))
      random_model_str = ""
      for (r_t in random_terms){random_model_str = paste0(random_model_str, " + (1|", r_t, ")")}  
      fml_fit = paste0(main_model_str, sv_model_str, random_model_str)
      fml_fixed = paste0(main_model_str, sv_model_str)

      if (distro_family == "poisson"){
        regression_caller = glmer
        #control_args = glmerControl()
      } else if (distro_family == "NB"){
        regression_caller = glmer.nb
        #control_args = NULL
        sp = 1
      }else if (distro_family == "gaussian"){
        regression_caller = lmer
        #control_args = lmerControl()
        sp = 1
      }
    }

    message("Fitting ", method, " models with ", distro_family, " distribution: ", formula(fml_fit), " ...")
    models0 = future_apply(irm, 1, fit_single_model, meta_data = irc, approach = sp, fit_formula = fml_fit, offset_str = offset_str, regression_caller = regression_caller, distro = distro, defaultlink = defaultlink)
    
    message("Checking model status...")
    names(models0) = rownames(irm)
    model_failed = sapply(models0, is.null)
    if (method == "glmm"){
      model_warned = sapply(models0, FUN = function(x){
        if (is.null(x)){
          warned = FALSE
        } else {
          warned = (length(x@optinfo$conv$lme4$messages) > 0)
        }
        return(warned)
        }
      )
    } else {
      model_warned = sapply(models0, FUN = function(x){
        if (is.null(x)){
          warned = FALSE
        } else {
          warned = (x$converged == FALSE)
        }
        return(warned)
        }
      )
    }
    models = models0[which((!model_failed) & (!model_warned))]
    model_status = data.frame(ok = rep(0, length(models0)), warned = rep(0, length(models0)), failed = rep(0, length(models0)))
    rownames(model_status) = names(models0)
    model_status[which(model_warned), "warned"] = 1
    model_status[which(model_failed), "failed"] = 1
    model_status[(model_status$warned == 0) & (model_status$failed == 0), "ok"] = 1
    message(length(models), " out of ", nrow(irm), " entries successfully modeled (", length(which(model_warned))," with warnings).")

    message("Extracting model parameters...")
    m = model.matrix(formula(fml_fixed), data = irc)
    m = as.data.frame(m)
    m_colnames = colnames(m)

    beta_all = t(sapply(models, FUN = model_coef_autocomplete, expected_coef_names = m_colnames, extract_column = 1))
    stderr_all = t(sapply(models, FUN = model_coef_autocomplete, expected_coef_names = m_colnames, extract_column = 2))
    p_all = t(sapply(models, FUN = model_coef_autocomplete, expected_coef_names = m_colnames, extract_column = 4))
    resid_all = t(sapply(models, resid, type = "response"))
    q_all = apply(as.matrix(p_all), 2, p.adjust, method = "BH")
    if (nrow(p_all) == 1){q_all = rbind(c(), q_all)}
    vcov_all = merge_vcov_from_model_list(models, expected_vcov_names = m_colnames)
    rownames(beta_all) = names(models)
    rownames(stderr_all) = names(models)
    rownames(p_all) = names(models)
    rownames(resid_all) = names(models)
    colnames(resid_all) = colnames(irm)
    rownames(q_all) = names(models)

    if ((length(c(remove_terms, sv_terms)) > 0) || (method != "glm")){
      message("Correcting against unwanted effects: ", paste(c(remove_terms, sv_terms), collapse = ", ")," ...")
      for (keyword in remove_terms){
        remove_cols = c(remove_cols, which(grepl(paste0("^", keyword), colnames(m))))
      }
      for (svword in sv_terms){
        remove_cols = c(remove_cols, which(colnames(m) == svword))
      }
      remove_cols = unique(remove_cols)
      #X = m[, remove_cols]
      #X = as.matrix(X)
      #beta_sub = beta_all[, remove_cols]
      #beta_sub = as.matrix(beta_sub)
      #cleany = y - t(as.matrix(X) %*% t(beta_sub))
      m[, remove_cols] = 0
      cleanY = reverselink(as.matrix(beta_all) %*% t(as.matrix(m))) + resid_all
      cleanY[which(cleanY < 0)] = 0
      cleany = defaultlink(cleanY)
    } else {
      if (! apply_sva){
        message("Getting rescaled data as final `y`...")
        if ((distro_family == "gaussian") || (distro_family == "Gaussian")){
          irm_norm = irm
        } else {
          irm_norm = t(apply(irm, 1, FUN = function(x){return(x / as.vector(irc[[offset_str]]))}))
        }
      }
      cleany = as.matrix(defaultlink(irm_norm[names(models), , drop = FALSE]))
    }
    
    message("Collecting results...")
    rownames(cleany) = names(models)
    args_all = list(method = method,
                    estimate_offset = estimate_offset,
                    offset_str = offset_str,
                    apply_sva = apply_sva,
                    n_sv = n_sv,
                    test_terms = test_terms,
                    remove_terms = remove_terms,
                    random_terms = random_terms,
                    distro_family = distro_family)
    
    if (save_models){
      glm_res = new("GLMres", y = irm, colData = irc, design = formula(fml_fit), beta = as.matrix(beta_all), stderr = as.matrix(stderr_all), pval = as.matrix(p_all), fdr = as.matrix(q_all), vcov = vcov_all, resid = as.matrix(resid_all), deSV = as.matrix(cleany), status = as.matrix(model_status), args = args_all, models = models0)
    } else {
      glm_res = new("GLMres", y = irm, colData = irc, design = formula(fml_fit), beta = as.matrix(beta_all), stderr = as.matrix(stderr_all), pval = as.matrix(p_all), fdr = as.matrix(q_all), vcov = vcov_all, resid = as.matrix(resid_all), deSV = as.matrix(cleany), status = as.matrix(model_status), args = args_all, models = list())
    }
    
    message("Done.")
    return(glm_res)
  }
}
