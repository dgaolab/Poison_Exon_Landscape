#!/usr/bin/env python3

import argparse
import sys, os
from subprocess import Popen, PIPE
from Bio.Seq import Seq
import pyBigWig
import numpy as np


def get_ir_coords(pid):
    #chr1/1084384/1084480/C1orf159/ENSG00000131591/-
    ch, start, end, gn, gid, st = pid.split("/")
    key_tmp = f'{ch}_{start}_{end}_{st}'
    
    return ('NA', 'NA', key_tmp)


def get_exonSkip_coords(pid):
    #ENSG00000160087@UBE2J2@chr1:-:1263387^1266097&1266291^1267861|1263387^1267861
    gid, gn, tmp = pid.split("@")
    ch, st, coords = tmp.split(":")
    t1 = coords.split("&")[0]
    t2 = coords.split("&")[1].split("|")[0]
    t3 = coords.split("|")[1]
    r1 = f'{ch}_{t1.split("^")[0]}_{t1.split("^")[1]}_{st}'
    r2 = f'{ch}_{t2.split("^")[0]}_{t2.split("^")[1]}_{st}'
    r3 = f'{ch}_{t3.split("^")[0]}_{t3.split("^")[1]}_{st}'

    return (r1, r2, r3)


def get_altExon_coords(pid):
    #ENSG00000188976@NOC2L@chr1:-:953893^954003&954524^955922|953893^954003&954083^955922
    gid, gn, tmp = pid.split("@")
    ch, st, coords = tmp.split(":")
    t1 = coords.split("&")[0]
    t2 = coords.split("&")[1].split("|")[0]
    t3 = coords.split("|")[1].split("&")[0]
    t4 = coords.split("&")[-1]

    if t1 == t3:
        a1 = t2
        a2 = t4
    elif t2 == t4:
        a1 = t1
        a2 = t3
    else:
        sys.exit("Error: unrealistic alternative-length exon.")

    r1 = f'{ch}_{a1.split("^")[0]}_{a1.split("^")[1]}_{st}'
    r3 = f'{ch}_{a2.split("^")[0]}_{a2.split("^")[1]}_{st}'

    return (r1, "NA", r3)


def get_seq(fa, ch, start, end, st):
    if "chr" in ch == False:
        ch_new = "chr" + ch
    else:
        ch_new = ch
    if (start != 0) and (end != 0):
        cmd = f'samtools faidx {fa} {ch_new}:{start}-{end}'  #1-based coordinates
        prc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out,stat = prc.communicate()
        seq0 = "".join(out.decode().split("\n")[1:])
        if st == "-":
            seq0 = str(Seq(seq0).reverse_complement())
    else:
        seq0 = ""

    return(seq0)


def get_attribute(info_str, attribute_key):
    if attribute_key in info_str:
        return info_str.split((attribute_key + " "))[1].split(";")[0][1:-1]
    else:
        sys.exit(f"'{attribute_key}' was not found in the attributes.")


def get_phylop(pybw, ch, start, end):
    avg_p = pybw.stats(ch, start, end, type = "mean")[0]
    if avg_p == None:
        avg_p = 0
    values = pybw.values(ch, start, end)
    values = np.array(values)
    valid = values[~np.isnan(values)]
    if len(valid) > 0:
        med_p = np.median(valid)
    else:
        med_p = 0

    return (avg_p, med_p)


def make_triplet_es_mode(i1, i2, ch, st, gene_id, gene_name):
    intron1_start = i1[0] + 1
    intron1_end = i1[1]
    intron2_start = i2[0] + 1
    intron2_end = i2[1]
    #ENSG00000000003@TSPAN6@chrX:-:100629987^100630758&100630867^100632484|100629987^100632484
    splice_id = f'{gene_id}@{gene_name}@{ch}:{st}:{intron1_start}^{intron1_end}&{intron2_start}^{intron2_end}|{intron1_start}^{intron2_end}'

    return splice_id


def make_triplet_ea_mode(intron1_l, intron2_l, intron1_s, intron2_s, ch, st, gene_id, gene_name):
    le_intron1_start = intron1_l[0] + 1
    le_intron1_end   = intron1_l[1]
    le_intron2_start = intron2_l[0] + 1
    le_intron2_end   = intron2_l[1]
    se_intron1_start = intron1_s[0] + 1
    se_intron1_end   = intron1_s[1]
    se_intron2_start = intron2_s[0] + 1
    se_intron2_end   = intron2_s[1]
    #ENSG00000000003@TSPAN6@chrX:-:100629987^100630758&100630867^100632484|100630758&100630867
    splice_id = f'{gene_id}@{gene_name}@{ch}:{st}:{le_intron1_start}^{le_intron1_end}&{le_intron2_start}^{le_intron2_end}|{se_intron1_start}^{se_intron1_end}&{se_intron2_start}^{se_intron2_end}'

    return splice_id


def contain_continuous_introns(i1, i2, tx, tx_dic):
    hit = 0
    if (i1 in tx_dic[tx]["introns"]) and (i2 in tx_dic[tx]["introns"]):
        i1_idx = tx_dic[tx]["introns"].index(i1)
        i2_idx = tx_dic[tx]["introns"].index(i2)
        if (i1_idx + 1) == i2_idx:
            hit = 1

    return hit


def locate_upstream_exon(intron, tx_cd, tx_dic, st):
    idx = tx_dic[tx_cd]["introns"].index(intron)
    if st == "+":
        up_ex = tx_dic[tx_cd]["struct"].split("|")[idx]
    else:
        up_ex = tx_dic[tx_cd]["struct"].split("|")[idx + 1]

    return up_ex


def locate_middle_exon(i1, i2, tx_cd, tx_dic):
    i2_idx = tx_dic[tx_cd]["introns"].index(i2)
    mid_ex = tx_dic[tx_cd]["struct"].split("|")[i2_idx]

    return mid_ex


def tandem_exons(gene_lines):
    tx_dic = {}
    exon_dic = {}
    intron_pairs = []

    line_idx = 0
    for x in gene_lines:
        ch, src, tp, start, end, null1, st, null2, info = x.rstrip().split("\t")
        if tp == "gene":
            gn_id = get_attribute(info, "gene_id").split(".")[0]
            if "gene_name" in info:
                gn_nm = get_attribute(info, "gene_name")
            else:
                gn_nm = "Novel_Gene"
            tx_current = ""
        
        if tp == "transcript":
            if tx_current != "":
                if st == "-":
                    tx_exons.reverse()      # DNA-level 5' -> 3' order, not RNA
                    exon_coords.reverse()   # DNA-level 5' -> 3' order, not RNA
                    cds_coords.reverse()    # DNA-level 5' -> 3' order, not RNA
                    cds_phases.reverse()    # DNA-level 5' -> 3' order, not RNA
                tx_struct = "|".join(tx_exons)
                intron_coords = [(exon_coords[i][1], exon_coords[i + 1][0]) for i in range(len(exon_coords) - 1)]
                if len(exon_coords) >= 2:
                    if st == "+":
                        last1exon_len = abs(exon_coords[-1][0] - exon_coords[-1][1])
                        last2exon_len = last1exon_len + abs(exon_coords[-2][0] - exon_coords[-2][1])
                    else:
                        last1exon_len = abs(exon_coords[0][0] - exon_coords[0][1])
                        last2exon_len = last1exon_len + abs(exon_coords[1][0] - exon_coords[1][1])
                else:
                    last1exon_len = tx_len
                    last2exon_len = tx_len
                for i in range(len(intron_coords) - 1):
                    intron_set = (intron_coords[i][0], intron_coords[i][1], intron_coords[i + 1][0], intron_coords[i + 1][1])
                    if (intron_set in intron_pairs) == False:
                        intron_pairs.append(intron_set)
                tx_phase_list = [p for p, (s, e) in zip(cds_phases, cds_coords) if (s, e) != (0, 0)]
                if st == "+":
                    tx_phase = tx_phase_list[0] if tx_phase_list else 0
                else:
                    tx_phase = tx_phase_list[-1] if tx_phase_list else 0
                tx_dic[(tx_current + "|" + tx_tp)] = {"length": tx_len, "struct":tx_struct, "strand": st, "CDS":cds_coords, "CDS_phases":cds_phases, "introns":intron_coords, "TSS": tss, "TES": tes, "last2exon_length": last2exon_len, "last1exon_length": last1exon_len, "UTR5": utr5, "UTR5_length": utr5_len, "UTR3": utr3, "UTR3_length": utr3_len, "Start_Codon": start_codon, "Stop_Codon": stop_codon, "cds_start_NF": cds_start_nf, "cds_end_NF": cds_end_nf, "mRNA_start_NF": mRNA_start_nf, "mRNA_end_NF": mRNA_end_nf, "ORF_shift": tx_phase}

            tx_id = get_attribute(info, "transcript_id")
            tx_current = tx_id
            if "transcript_biotype" in info:
                tx_tp = get_attribute(info, "transcript_biotype")
            else:
                tx_tp = get_attribute(info, "transcript_type")
            if st == "+":
                tss = int(start) - 1    # 0-based coordinates
                tes = int(end) - 1      # 0-based coordinates
            else:
                tss = int(end) - 1      # 0-based coordinates
                tes = int(start) - 1    # 0-based coordinates
            cds_start_nf = "cds_start_NF" in info
            cds_end_nf = "cds_end_NF" in info
            mRNA_start_nf = "mRNA_start_NF" in info
            mRNA_end_nf = "mRNA_end_NF" in info
            tx_exons = []
            exon_coords = []            # 0-based coordinates
            cds_coords = []             # 0-based coordinates
            cds_phases = []             # CDS frame/phase from GTF column 8
            start_codon = (-1, -1)      # 0-based coordinates
            stop_codon = (-1, -1)       # 0-based coordinates
            utr5 = []                   # 0-based coordinates
            utr3 = []                   # 0-based coordinates
            tx_len = 0
            utr5_len = 0
            utr3_len = 0

        if tp == "exon":
            exon_tx_id = get_attribute(info, "transcript_id")
            if exon_tx_id != tx_current:
                sys.exit("'transcript_id' in 'transcript' and 'exon' lines doesn't match.")
            exon_coords.append((int(start) - 1, int(end)))            # convert exon coordinates to 0-based
            exon_id = get_attribute(info, "exon_id")
            if (exon_id in exon_dic) == False:
                exon_dic[exon_id] = (int(start) - 1, int(end))        # convert exon coordinates to 0-based
            exon_len = int(end) - int(start) + 1
            tx_len += exon_len
            tx_exons.append(exon_id)
            if (line_idx + 1) < len(gene_lines):
                tmp = gene_lines[line_idx + 1].rstrip().split("\t")
                if tmp[2] == "CDS":
                    cds_coords.append((int(tmp[3]) - 1, int(tmp[4]))) # convert CDS coordinates to 0-based
                    cds_phases.append(int(tmp[7]) if tmp[7] != "." else 0)
                else:
                    cds_coords.append((0, 0))
                    cds_phases.append(0)
            else:
                cds_coords.append((0, 0))
                cds_phases.append(0)
        
        if tp in ["start_codon", "stop_codon"]:
            current_tx_id = get_attribute(info, "transcript_id")
            if current_tx_id != tx_current:
                sys.exit(f"'transcript_id' in 'transcript' and start/stop_codon lines doesn't match.")
            if tp == "start_codon":
                start_codon = (int(start) - 1, int(end))
            elif tp == "stop_codon":
                stop_codon = (int(start) - 1, int(end))

        if tp in ["UTR", "five_prime_utr", "three_prime_utr"]:
            current_tx_id = get_attribute(info, "transcript_id")
            if current_tx_id != tx_current:
                sys.exit(f"'transcript_id' in 'transcript' and UTR lines doesn't match.")
            cds_start_list = [x1 for x1, x2 in cds_coords if x1 > 0]
            if len(cds_start_list) > 0:
                cds_start = min(cds_start_list)
            else:
                cds_start = 0
            cds_end_list = [x2 for x1, x2 in cds_coords if x2 > 0]
            if len(cds_end_list) > 0:
                cds_end = max(cds_end_list)
            else:
                cds_end = 0
            if tp == "five_prime_utr":
                utr5.append((int(start) - 1, int(end)))
                utr5_len += (int(end) - int(start) + 1)
            elif tp == "three_prime_utr":
                utr3.append((int(start) - 1, int(end)))
                utr3_len += (int(end) - int(start) + 1)
            elif tp == "UTR":
                if st == "+":
                    if (int(end) <= cds_start) and (cds_start > 0):
                        utr5.append((int(start) - 1, int(end)))
                        utr5_len += (int(end) - int(start) + 1)
                    elif (int(start) >= cds_end) and (cds_end > 0):
                        utr3.append((int(start) - 1, int(end)))
                        utr3_len += (int(end) - int(start) + 1)
                    else:
                        print(cds_start, cds_end, start, end)
                        sys.exit("Error: Ambiguous UTR.")
                else:
                    if (int(end) <= cds_start) and (cds_start > 0):
                        utr3.append((int(start) - 1, int(end)))
                        utr3_len += (int(end) - int(start) + 1)
                    elif (int(start) >= cds_end) and (cds_end > 0):
                        utr5.append((int(start) - 1, int(end)))
                        utr5_len += (int(end) - int(start) + 1)
                    else:
                        sys.exit("Error: Ambiguous UTR.")
        
        line_idx += 1

    #record the last transcript in a gene
    if st == "-":
        tx_exons.reverse()       # DNA-level 5' -> 3' order, not RNA
        exon_coords.reverse()    # DNA-level 5' -> 3' order, not RNA
        cds_coords.reverse()     # DNA-level 5' -> 3' order, not RNA
        cds_phases.reverse()     # DNA-level 5' -> 3' order, not RNA
    tx_struct = "|".join(tx_exons)
    intron_coords = [(exon_coords[i][1], exon_coords[i + 1][0]) for i in range(len(exon_coords) - 1)]
    if len(exon_coords) >= 2:
        if st == "+":
            last1exon_len = abs(exon_coords[-1][0] - exon_coords[-1][1])
            last2exon_len = last1exon_len + abs(exon_coords[-2][0] - exon_coords[-2][1])
        else:
            last1exon_len = abs(exon_coords[0][0] - exon_coords[0][1])
            last2exon_len = last1exon_len + abs(exon_coords[1][0] - exon_coords[1][1])
    else:
        last1exon_len = tx_len
        last2exon_len = tx_len
    for i in range(len(intron_coords) - 1):
        intron_set = (intron_coords[i][0], intron_coords[i][1], intron_coords[i + 1][0], intron_coords[i + 1][1])
        if (intron_set in intron_pairs) == False:
            intron_pairs.append(intron_set)
    tx_phase_list = [p for p, (s, e) in zip(cds_phases, cds_coords) if (s, e) != (0, 0)]
    if st == "+":
        tx_phase = tx_phase_list[0] if tx_phase_list else 0
    else:
        tx_phase = tx_phase_list[-1] if tx_phase_list else 0
    tx_dic[(tx_current + "|" + tx_tp)] = {"length": tx_len, "struct": tx_struct, "strand": st, "CDS": cds_coords, "CDS_phases": cds_phases, "introns": intron_coords, "TSS": tss, "TES": tes, "last2exon_length": last2exon_len, "last1exon_length": last1exon_len, "UTR5": utr5, "UTR5_length": utr5_len, "UTR3": utr3, "UTR3_length": utr3_len, "Start_Codon": start_codon, "Stop_Codon": stop_codon, "cds_start_NF": cds_start_nf, "cds_end_NF": cds_end_nf, "mRNA_start_NF": mRNA_start_nf, "mRNA_end_NF": mRNA_end_nf, "ORF_shift": tx_phase}

    return (ch, st, gn_id, gn_nm, intron_pairs, tx_dic, exon_dic)


def get_pe_utr_overlap(start, end, tx_id, tx_dic):
    if tx_dic[tx_id]["UTR5"] == []:
        hit_utr5 = 0
    else:
        utr5_start = min([x1 for x1, x2 in tx_dic[tx_id]["UTR5"]])
        utr5_end = max([x2 for x1, x2 in tx_dic[tx_id]["UTR5"]])
        hit_utr5 = int(max(start, utr5_start) < min(end, utr5_end))
    if tx_dic[tx_id]["UTR3"] == []:
        hit_utr3 = 0
    else:
        utr3_start = min([x1 for x1, x2 in tx_dic[tx_id]["UTR3"]])
        utr3_end = max([x2 for x1, x2 in tx_dic[tx_id]["UTR3"]])
        hit_utr3 = int(max(start, utr3_start) < min(end, utr3_end))
    hit_utr = max([hit_utr5, hit_utr3])

    return hit_utr


def if_ptc_by_ins(fa, ch, st, tx_cd, tx_dic, ins_start, ins_end, upstream_exon):
    stop_codon = ["TAG", "TGA", "TAA"]
    
    ins_seq = get_seq(fa, ch, ins_start + 1, ins_end, st)   #convert to 1-based coordinates
    seg_seq_list = [get_seq(fa, ch, pos1 + 1, pos2, st) for pos1, pos2 in tx_dic[tx_cd]["CDS"]] #convert to 1-based coordinates
    seg_names = tx_dic[tx_cd]["struct"].split("|")
    if st == "-":
        seg_seq_list.reverse()
        seg_names.reverse()

    upstream_exon_idx = seg_names.index(upstream_exon)
    seg_seq_list_new = seg_seq_list[:(upstream_exon_idx + 1)] + [ins_seq] + seg_seq_list[(upstream_exon_idx + 1):]
    seg_seq_new = "".join(seg_seq_list_new)
    
    phase = tx_dic[tx_cd]["ORF_shift"]
    orf_codons = [seg_seq_new[ii:(ii+3)] for ii in range(phase, len(seg_seq_new)) if ((ii - phase) % 3) == 0]
    orf_codons_stop = [1 if codon in stop_codon else 0 for codon in orf_codons]
    ptc_found = max(orf_codons_stop)
    
    if ptc_found == 1:
        ptc_found_index = orf_codons_stop.index(ptc_found)
        ptc_pos = phase + ptc_found_index * 3
    else:
        ptc_pos = -1

    return ptc_pos


def if_ptc_by_del(fa, ch, st, tx_cd, tx_dic, short_start, short_end, target_exon):
    stop_codon = ["TAG", "TGA", "TAA"]
    
    seg_names = tx_dic[tx_cd]["struct"].split("|")
    target_exon_idx = seg_names.index(target_exon)
    target_exon_start = tx_dic[tx_cd]["CDS"][target_exon_idx][0]  #0-based coordinate
    target_exon_end   = tx_dic[tx_cd]["CDS"][target_exon_idx][1]  #0-based coordinate

    target_exon_cd_len = target_exon_end - target_exon_start
    del_cd_len = short_end - short_start

    if del_cd_len >= target_exon_cd_len:
        ptc_pos = -1
    else:
        replace_seq = get_seq(fa, ch, short_start + 1, short_end, st)   #convert to 1-based coordinates
        seg_seq_list = [get_seq(fa, ch, pos1 + 1, pos2, st) for pos1, pos2 in tx_dic[tx_cd]["CDS"]]  #convert to 1-based coordinates
        if st == "-":
            seg_seq_list.reverse()
            seg_names.reverse()
            target_exon_idx = len(seg_names) - target_exon_idx - 1

        seg_seq_list_new = seg_seq_list[:target_exon_idx] + [replace_seq] + seg_seq_list[(target_exon_idx + 1):]
        seg_seq_new = "".join(seg_seq_list_new)
        
        phase = tx_dic[tx_cd]["ORF_shift"]
        orf_codons = [seg_seq_new[ii:(ii+3)] for ii in range(phase, len(seg_seq_new)) if ((ii - phase) % 3) == 0]
        orf_codons_stop = [1 if codon in stop_codon else 0 for codon in orf_codons]
        ptc_found = max(orf_codons_stop)
        
        if ptc_found == 1:
            ptc_found_index = orf_codons_stop.index(ptc_found)
            ptc_pos = phase + ptc_found_index * 3
        else:
            ptc_pos = -1

    return ptc_pos


def get_ptc_relative_position(pe_len, ptc_pos, tx_id, tx_dic):
    # Determine phase offset of first CDS in 5'->3' order
    phase = tx_dic[tx_id]["ORF_shift"]

    # Calculate PTC's distance to TSS
    # `ptc_pos` is relative to the start of CDS
    # Note: when mRNA_start_NF or cds_start_NF, TSS is unknown, we still calculate this number based on the known ORF start 
    dis_tss = ptc_pos + tx_dic[tx_id]["UTR5_length"] + phase

    # We need to convert all genomic coordinates (e.g. exon, UTR) to transcriptomic coordinates relative to TSS first
    # Note: when mRNA_end_NF or cds_end_NF, TES is unknown, we still calculate this number based on the known end
    ref_pos_last2exon = tx_dic[tx_id]["length"] + pe_len - tx_dic[tx_id]["last2exon_length"]
    ref_pos_last1exon = tx_dic[tx_id]["length"] + pe_len - tx_dic[tx_id]["last1exon_length"]

    # For last-exon checks, we need a valid dis_tss; use computed value even if flagged -1 for TSS
    dis_tss_for_check = ptc_pos + tx_dic[tx_id]["UTR5_length"]

    # whether PTC is downstream to the start of the 2nd last exon
    if (dis_tss_for_check + 3) < ref_pos_last2exon:
        hit_2nd_last = 0
    else:
        hit_2nd_last = 1

    # whether PTC is downstream to the start of the last exon
    if (dis_tss_for_check + 3) < ref_pos_last1exon:
        hit_last = 0
    else:
        hit_last = 1

    return (dis_tss, hit_2nd_last, hit_last)


def get_pe_ptc_overlap(upstream_exon, pe_len, ptc_pos, tx_id, tx_dic):
    seg_len_list = [(pos2 - pos1) for pos1, pos2 in tx_dic[tx_id]["CDS"]] #convert to 1-based coordinates
    seg_names = tx_dic[tx_id]["struct"].split("|")
    if tx_dic[tx_id]["strand"] == "-":
        seg_len_list.reverse()
        seg_names.reverse()

    upstream_exon_idx = seg_names.index(upstream_exon)
    pe_start_pos = sum(seg_len_list[:(upstream_exon_idx + 1)])
    pe_end_pos = pe_start_pos + pe_len
    if ptc_pos <= pe_end_pos:
        hit_pe = 1
    else:
        hit_pe = 0

    return hit_pe


def judge_exon_skipping(fa, ch, st, gn_id, gn_nm, intron_pairs, tx_dic, writer0, writer1, pybw):
    for intron_pair in intron_pairs:
        i1 = intron_pair[0:2]
        i2 = intron_pair[2:]
        i0 = (i1[0], i2[1])

        incl_tx = [tx for tx in tx_dic if contain_continuous_introns(i1, i2, tx, tx_dic)]
        excl_tx = [tx for tx in tx_dic if i0 in tx_dic[tx]["introns"]]

        both_forms = 0
        only_inclusion_form = 0
        if (len(incl_tx) > 0) and (len(excl_tx) > 0):  #both exon-included and exon-excluded transcripts exist
            both_forms = 1
        if (len(incl_tx) > 0) and (len(excl_tx) == 0):
            only_inclusion_form = 1
        if (len(incl_tx) == 0) and (len(excl_tx) > 0):
            sys.exit("Error: exons always excluded from a gene's transcripts shouldn't exist in the annotation.")

        if both_forms:
            splice_id = make_triplet_es_mode(i1, i2, ch, st, gn_id, gn_nm)
            region1, region2, region3 = get_exonSkip_coords(splice_id)
            tx_with_str = ";".join(incl_tx)
            tx_without_str = ";".join(excl_tx)
            writer0.write(f'{splice_id}\tref_exon_skipping\t{tx_with_str}\t{tx_without_str}\t{region1}\t{region2}\t{region3}\n')

            e2_start = i1[1]               #middle exon start; 0-based coordinate
            e2_end   = i2[0]               #middle exon end;   0-based coordinate
            exon_len = e2_end - e2_start 
            
            if writer1 != None:
                #Seems we have to find:
                #1) 'protein_coding' exon-EXCLUDED transcript exists in the annotation;
                #2) exon-included transcript exists in the annotation and NONE of them should be 'protein_coding';
                #3) exon length less than a threshold (default: 300bp)
                #4) Insert exon into the 'protein_coding' exon-EXCLUDED transcript will lead to PTC in the new transcript;
                if ("protein_coding" in tx_without_str) and (("protein_coding" in tx_with_str) == False) and ((exon_len) < 300):
                    doublet_coding = [xx for xx in excl_tx if "protein_coding" in xx]
                    ptc_res = []
                    for tx_cd in doublet_coding:
                        pe_utr = get_pe_utr_overlap(e2_start, e2_end, tx_cd, tx_dic)
                        if pe_utr == 0:    # Ensure PE is not inserted into UTR
                            upstream_exon = locate_upstream_exon(i0, tx_cd, tx_dic, st)
                            ptc_hit = if_ptc_by_ins(fa, ch, st, tx_cd, tx_dic, e2_start, e2_end, upstream_exon)
                            if ptc_hit > 0:
                                ptc_tss, ptc_2nd_last, ptc_last = get_ptc_relative_position(exon_len, ptc_hit, tx_cd, tx_dic)
                                ptc_pe = get_pe_ptc_overlap(upstream_exon, exon_len, ptc_hit, tx_cd, tx_dic)
                                ptc_res.append((ptc_hit, tx_cd, str(ptc_tss), str(ptc_2nd_last), str(ptc_last), str(ptc_pe), str(int(tx_dic[tx_cd]["mRNA_start_NF"])), str(int(tx_dic[tx_cd]["cds_start_NF"])), str(int(tx_dic[tx_cd]["mRNA_end_NF"])), str(int(tx_dic[tx_cd]["cds_end_NF"])), str(tx_dic[tx_cd]["ORF_shift"])))

                    if ptc_res != []:
                        ptc_tx = [trans for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                        ptc_tx_str = ";".join(ptc_tx)
                        d2tss = [dis_tss for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                        d2tss_str = ";".join(d2tss)
                        hit_tail = [on_tail for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                        hit_tail_str = ";".join(hit_tail)
                        hit_last = [on_last for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                        hit_last_str = ";".join(hit_last)
                        hit_pe = [on_pe for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                        hit_pe_str = ";".join(hit_pe)
                        ms_nf_str = ";".join([ms_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                        cs_nf_str = ";".join([cs_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                        me_nf_str = ";".join([me_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                        ce_nf_str = ";".join([ce_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                        shift_str = ";".join([shift for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                        if pybw != None:
                            pe_avg_phylop, pe_med_phylop = get_phylop(pybw, ch, e2_start, e2_end)
                        else:
                            pe_avg_phylop = -999
                            pe_med_phylop = -999
                        writer1.write(f'{splice_id}\tref_exon_poison_by_inclusion\t{exon_len}\t{pe_avg_phylop}\t{pe_med_phylop}\t{ptc_tx_str}\t{shift_str}\t{ms_nf_str}\t{cs_nf_str}\t{me_nf_str}\t{ce_nf_str}\t{d2tss_str}\t{hit_tail_str}\t{hit_last_str}\t{hit_pe_str}\t{region1}\t{region2}\t{region3}\n')

        if only_inclusion_form: # no alternative form, always included. 
            splice_id = make_triplet_es_mode(i1, i2, ch, st, gn_id, gn_nm)
            region1, region2, region3 = get_exonSkip_coords(splice_id)
            tx_with_str = ";".join(incl_tx)
            writer0.write(f'{splice_id}\tref_exon_alwaysIn\t{tx_with_str}\tNA\t{region1}\t{region2}\t{region3}\n')

    return


def judge_exon_altLen(fa, ch, st, gn_id, gn_nm, intron_pairs, tx_dic, writer0, writer1, pybw):
    w0 = []
    w1 = []
    for intron_pair1 in intron_pairs:
        ir1 = intron_pair1[0:2]
        ir2 = intron_pair1[2:]
        er2 = (ir1[1], ir2[0])        #middle exon (start, end); 0-based coordinates
        for intron_pair2 in intron_pairs:
            both_forms = 0
            triplet_type = ""

            it1 = intron_pair2[0:2]
            it2 = intron_pair2[2:]
            et2 = (it1[1], it2[0])    #middle exon (start, end); 0-based coordinates

            if (intron_pair1 != intron_pair2) and (ir1[0] == it1[0]) and (ir2[1] == it2[1]):
                if (er2[0] == et2[0]) and (er2[1] != et2[1]):
                    both_forms = 1
                    if st == "+":
                        triplet_type = "ref_exon_alt3"
                    elif st == "-":
                        triplet_type = "ref_exon_alt5"
                    if er2[1] > et2[1]:
                        il1, el2, il2 = (ir1, er2, ir2)
                        is1, es2, is2 = (it1, et2, it2)
                    else:
                        il1, el2, il2 = (it1, et2, it2)
                        is1, es2, is2 = (ir1, er2, ir2)
                    diff_start = es2[1]
                    diff_end   = el2[1]

                if (er2[0] != et2[0]) and (er2[1] == et2[1]):
                    both_forms = 1
                    if st == "+":
                        triplet_type = "ref_exon_alt5"
                    elif st == "-":
                        triplet_type = "ref_exon_alt3"
                    if er2[0] < et2[0]:
                        il1, el2, il2 = (ir1, er2, ir2)
                        is1, es2, is2 = (it1, et2, it2)
                    else:
                        il1, el2, il2 = (it1, et2, it2)
                        is1, es2, is2 = (ir1, er2, ir2)
                    diff_start = el2[0]          #0-based coordinate
                    diff_end   = es2[0]          #0-based coordinate

            if both_forms:
                splice_id = make_triplet_ea_mode(il1, il2, is1, is2, ch, st, gn_id, gn_nm)
                region1, region2, region3 = get_altExon_coords(splice_id)
                le_tx = [tx for tx in tx_dic if contain_continuous_introns(il1, il2, tx, tx_dic)]
                se_tx = [tx for tx in tx_dic if contain_continuous_introns(is1, is2, tx, tx_dic)]
                tx_le_str = ";".join(le_tx)
                tx_se_str = ";".join(se_tx)
                w0_line = f'{splice_id}\t{triplet_type}\t{tx_le_str}\t{tx_se_str}\t{region1}\t{region2}\t{region3}\n'
                if (w0_line in w0) == False:
                    writer0.write(w0_line)
                    w0.append(w0_line)

                if writer1 != None:
                    exon_l_len = el2[1] - el2[0]  #0-based coordinate
                    exon_s_len = es2[1] - es2[0]  #0-based coordinate
                    if ("protein_coding" in tx_le_str) and (("protein_coding" in tx_se_str) == False):
                        coding_switchable = 1
                    elif ("protein_coding" in tx_se_str) and (("protein_coding" in tx_le_str) == False):
                        coding_switchable = 1
                    else:
                        coding_switchable = 0
                    
                    if (coding_switchable == 1) and (min(exon_l_len, exon_s_len) < 300):
                        #Either the short or the long form of the exon must contain an in-frame PTC
                        #The other form has to be in a protein_coding transcript
                        tx_le_cd = [xx for xx in le_tx if "protein_coding" in xx]
                        tx_se_cd = [xx for xx in se_tx if "protein_coding" in xx]

                        ptc_res = []
                        for tx_cd in tx_se_cd:
                            pe_utr = get_pe_utr_overlap(el2[0], el2[1], tx_cd, tx_dic)
                            if pe_utr == 0:    # Ensure PE is not inserted into UTR
                                if st == "+":
                                    if triplet_type == "ref_exon_alt5":
                                        upstream_exon = locate_upstream_exon(is1, tx_cd, tx_dic, st)
                                    elif triplet_type == "ref_exon_alt3":
                                        upstream_exon = locate_middle_exon(is1, is2, tx_cd, tx_dic)
                                else:
                                    if triplet_type == "ref_exon_alt5":
                                        upstream_exon = locate_upstream_exon(is2, tx_cd, tx_dic, st)
                                    elif triplet_type == "ref_exon_alt3":
                                        upstream_exon = locate_middle_exon(is1, is2, tx_cd, tx_dic)
                                ptc_hit = if_ptc_by_ins(fa, ch, st, tx_cd, tx_dic, diff_start, diff_end, upstream_exon)
                                if (ptc_hit > 0) and (exon_l_len < 300):
                                    ptc_tss, ptc_2nd_last, ptc_last = get_ptc_relative_position(exon_l_len, ptc_hit, tx_cd, tx_dic)
                                    if st == "+":
                                        previous_exon = locate_upstream_exon(is1, tx_cd, tx_dic, st)
                                    else:
                                        previous_exon = locate_upstream_exon(is2, tx_cd, tx_dic, st)
                                    ptc_pe = get_pe_ptc_overlap(previous_exon, exon_l_len, ptc_hit, tx_cd, tx_dic)
                                    ptc_res.append((ptc_hit, tx_cd, str(ptc_tss), str(ptc_2nd_last), str(ptc_last), str(ptc_pe), str(int(tx_dic[tx_cd]["mRNA_start_NF"])), str(int(tx_dic[tx_cd]["cds_start_NF"])), str(int(tx_dic[tx_cd]["mRNA_end_NF"])), str(int(tx_dic[tx_cd]["cds_end_NF"])), str(tx_dic[tx_cd]["ORF_shift"])))

                        if ptc_res != []:
                            ptc_tx = [trans for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            ptc_tx_str = ";".join(ptc_tx)
                            d2tss = [dis_tss for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            d2tss_str = ";".join(d2tss)
                            hit_tail = [on_tail for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_tail_str = ";".join(hit_tail)
                            hit_last = [on_last for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_last_str = ";".join(hit_last)
                            hit_pe = [on_pe for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_pe_str = ";".join(hit_pe)
                            ms_nf_str = ";".join([ms_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            cs_nf_str = ";".join([cs_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            me_nf_str = ";".join([me_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            ce_nf_str = ";".join([ce_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            shift_str = ";".join([shift for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            if pybw != None:
                                pe_avg_phylop, pe_med_phylop = get_phylop(pybw, ch, el2[0], el2[1])
                            else:
                                pe_avg_phylop = -999
                                pe_med_phylop = -999
                            w1_line = f'{splice_id}\t{triplet_type + "_poison_by_long"}\t{exon_l_len}\t{pe_avg_phylop}\t{pe_med_phylop}\t{ptc_tx_str}\t{shift_str}\t{ms_nf_str}\t{cs_nf_str}\t{me_nf_str}\t{ce_nf_str}\t{d2tss_str}\t{hit_tail_str}\t{hit_last_str}\t{hit_pe_str}\t{region1}\t{region2}\t{region3}\n'
                            if (w1_line in w1) == False:
                                writer1.write(w1_line)
                                w1.append(w1_line)

                        ptc_res = []
                        for tx_cd in tx_le_cd:
                            pe_utr = get_pe_utr_overlap(es2[0], es2[1], tx_cd, tx_dic)
                            if pe_utr == 0:    # Ensure PE is not inserted into UTR
                                target_exon = locate_middle_exon(il1, il2, tx_cd, tx_dic)
                                short_start = es2[0]
                                short_end = es2[1]
                                ptc_hit = if_ptc_by_del(fa, ch, st, tx_cd, tx_dic, short_start, short_end, target_exon)
                                if (ptc_hit > 0) and (exon_s_len < 300):
                                    ptc_tss, ptc_2nd_last, ptc_last = get_ptc_relative_position(exon_s_len, ptc_hit, tx_cd, tx_dic)
                                    if st == "+":
                                        previous_exon = locate_upstream_exon(il1, tx_cd, tx_dic, st)
                                    else:
                                        previous_exon = locate_upstream_exon(il2, tx_cd, tx_dic, st)
                                    ptc_pe = get_pe_ptc_overlap(previous_exon, exon_s_len, ptc_hit, tx_cd, tx_dic)
                                    ptc_res.append((ptc_hit, tx_cd, str(ptc_tss), str(ptc_2nd_last), str(ptc_last), str(ptc_pe), str(int(tx_dic[tx_cd]["mRNA_start_NF"])), str(int(tx_dic[tx_cd]["cds_start_NF"])), str(int(tx_dic[tx_cd]["mRNA_end_NF"])), str(int(tx_dic[tx_cd]["cds_end_NF"])), str(tx_dic[tx_cd]["ORF_shift"])))

                        if ptc_res != []:
                            ptc_tx = [trans for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            ptc_tx_str = ";".join(ptc_tx)
                            d2tss = [dis_tss for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            d2tss_str = ";".join(d2tss)
                            hit_tail = [on_tail for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_tail_str = ";".join(hit_tail)
                            hit_last = [on_last for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_last_str = ";".join(hit_last)
                            hit_pe = [on_pe for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                            hit_pe_str = ";".join(hit_pe)
                            ms_nf_str = ";".join([ms_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            cs_nf_str = ";".join([cs_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            me_nf_str = ";".join([me_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            ce_nf_str = ";".join([ce_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            shift_str = ";".join([shift for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res])
                            if pybw != None:
                                pe_avg_phylop, pe_med_phylop = get_phylop(pybw, ch, es2[0], es2[1])
                            else:
                                pe_avg_phylop = -999
                                pe_med_phylop = -999
                            w1_line = f'{splice_id}\t{triplet_type + "_poison_by_short"}\t{exon_s_len}\t{pe_avg_phylop}\t{pe_med_phylop}\t{ptc_tx_str}\t{shift_str}\t{ms_nf_str}\t{cs_nf_str}\t{me_nf_str}\t{ce_nf_str}\t{d2tss_str}\t{hit_tail_str}\t{hit_last_str}\t{hit_pe_str}\t{region1}\t{region2}\t{region3}\n'
                            if (w1_line in w1) == False:
                                writer1.write(w1_line)
                                w1.append(w1_line)

    return
        

def judge_exon_ir(fa, ch, st, gn_id, gn_nm, exon_dic, tx_dic, writer0, writer1, pybw):
    w0 = []
    w1 = []
    ir_ok0 = {}
    ir_ok1 = {}
    for ex_0 in exon_dic:
        start0 = exon_dic[ex_0][0]
        end0   = exon_dic[ex_0][1]
        left_hit = []
        right_hit = []
        sub_pairs = []
        intron_sizes = []
        for ex_1 in exon_dic:
            start1 = exon_dic[ex_1][0]
            end1 = exon_dic[ex_1][1]
            if (start0 == start1) and (end1 < end0):
                left_hit.append(ex_1)
            if (start0 < start1) and (end1 == end0):
                right_hit.append(ex_1)
        if (len(left_hit) > 0) and (len(right_hit) > 0):
            for x in left_hit:
                left_end = exon_dic[x][1]
                for y in right_hit:
                    right_start = exon_dic[y][0]
                    if left_end < right_start:
                        dis = right_start - left_end
                        if dis >= 10:
                            sub_pairs.append((x, y))
                            intron_sizes.append(dis)

        if len(sub_pairs) > 0:
            for ep, sz in zip(sub_pairs, intron_sizes):
                both_forms = 0
                exon_left, exon_right = ep
    
                doublet = exon_left + "|" + exon_right
                ir_tx = [tx for tx in tx_dic if ex_0 in tx_dic[tx]["struct"]]
                ok_tx = [tx for tx in tx_dic if doublet in tx_dic[tx]["struct"]]
                
                if (len(ir_tx) > 0) and (len(ok_tx) > 0):
                    both_forms = 1

                if both_forms:
                    intron_start = exon_dic[exon_left][1]    # 0-based coordinate
                    intron_end   = exon_dic[exon_right][0]   # 0-based coordinate    
                    splice_id = f'{ch}/{intron_start + 1}/{intron_end}/{gn_nm}/{gn_id}/{st}' # convert to 1-based coordinates
                    region1, region2, region3 = get_ir_coords(splice_id)
            
                    tx_ir_str = ";".join(ir_tx)
                    tx_ok_str = ";".join(ok_tx)
                    if splice_id in ir_ok0:
                        for irt in ir_tx:
                            if (irt in ir_ok0[splice_id]["tx_ir_str"]) == False:
                                ir_ok0[splice_id]["tx_ir_str"] += (";" + irt)
                        for okt in ok_tx:
                            if (okt in ir_ok0[splice_id]["tx_ok_str"]) == False:
                                ir_ok0[splice_id]["tx_ok_str"] += (";" + okt)
                    else:
                        ir_ok0[splice_id] = {"tx_ir_str": tx_ir_str, "tx_ok_str": tx_ok_str, "region1": region1, "region2": region2, "region3": region3}

                    if writer1 != None:
                        if ("protein_coding" in tx_ok_str) and (("protein_coding" in tx_ir_str) == False) and (sz < 300):
                            tx_ok_cd = [xx for xx in ok_tx if "protein_coding" in xx]
                            ptc_res = []
                            for tx_cd in tx_ok_cd:
                                pe_utr = get_pe_utr_overlap(intron_start, intron_end, tx_cd, tx_dic)
                                if pe_utr == 0:    # Ensure PE is not inserted into UTR
                                    if st == "+":
                                        upstream_exon = exon_left
                                    else:
                                        upstream_exon = exon_right
                                    ptc_hit = if_ptc_by_ins(fa, ch, st, tx_cd, tx_dic, intron_start, intron_end, upstream_exon)
                                    if ptc_hit > 0:
                                        ptc_tss, ptc_2nd_last, ptc_last = get_ptc_relative_position(sz, ptc_hit, tx_cd, tx_dic)
                                        ptc_pe = get_pe_ptc_overlap(upstream_exon, sz, ptc_hit, tx_cd, tx_dic)
                                        ptc_res.append((ptc_hit, tx_cd, str(ptc_tss), str(ptc_2nd_last), str(ptc_last), str(ptc_pe), str(int(tx_dic[tx_cd]["mRNA_start_NF"])), str(int(tx_dic[tx_cd]["cds_start_NF"])), str(int(tx_dic[tx_cd]["mRNA_end_NF"])), str(int(tx_dic[tx_cd]["cds_end_NF"])), str(tx_dic[tx_cd]["ORF_shift"])))

                            if ptc_res != []:
                                ptc_tx = [trans for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                ptc_tx_str = ";".join(ptc_tx)
                                d2tss = [dis_tss for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                d2tss_str = ";".join(d2tss)
                                hit_tail = [on_tail for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                hit_tail_str = ";".join(hit_tail)
                                hit_last = [on_last for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                hit_last_str = ";".join(hit_last)
                                hit_pe = [on_pe for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                hit_pe_str = ";".join(hit_pe)
                                ms_nf_list = [ms_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                ms_nf_str = ";".join(ms_nf_list)
                                cs_nf_list = [cs_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                cs_nf_str = ";".join(cs_nf_list)
                                me_nf_list = [me_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                me_nf_str = ";".join(me_nf_list)
                                ce_nf_list = [ce_nf for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                ce_nf_str = ";".join(ce_nf_list)
                                shift_list = [shift for sc, trans, dis_tss, on_tail, on_last, on_pe, ms_nf, cs_nf, me_nf, ce_nf, shift in ptc_res]
                                shift_str = ";".join(shift_list)
                                if pybw != None:
                                    pe_avg_phylop, pe_med_phylop = get_phylop(pybw, ch, intron_start, intron_end)
                                else:
                                    pe_avg_phylop = -999
                                    pe_med_phylop = -999

                                if splice_id in ir_ok1:
                                    for ptct_idx in range(len(ptc_tx)):
                                        if (ptc_tx[ptct_idx] in ir_ok1[splice_id]["ptc_tx_str"]) == False:
                                            ir_ok1[splice_id]["ptc_tx_str"] += (";" + ptc_tx[ptct_idx])
                                            ir_ok1[splice_id]["d2tss_str"] += (";" + str(d2tss[ptct_idx]))
                                            ir_ok1[splice_id]["hit_tail_str"] += (";" + str(hit_tail[ptct_idx]))
                                            ir_ok1[splice_id]["hit_last_str"] += (";" + str(hit_last[ptct_idx]))
                                            ir_ok1[splice_id]["hit_pe_str"] += (";" + str(hit_pe[ptct_idx]))
                                            ir_ok1[splice_id]["ms_nf_str"] += (";" + str(ms_nf_list[ptct_idx]))
                                            ir_ok1[splice_id]["cs_nf_str"] += (";" + str(cs_nf_list[ptct_idx]))
                                            ir_ok1[splice_id]["me_nf_str"] += (";" + str(me_nf_list[ptct_idx]))
                                            ir_ok1[splice_id]["ce_nf_str"] += (";" + str(ce_nf_list[ptct_idx]))
                                            ir_ok1[splice_id]["shift_str"] += (";" + str(shift_list[ptct_idx]))
                                else:
                                    ir_ok1[splice_id] = {"pe_len": sz, "avg_phylop": pe_avg_phylop, "med_phylop": pe_med_phylop, "ptc_tx_str": ptc_tx_str, "d2tss_str": d2tss_str, "hit_tail_str": hit_tail_str, "hit_last_str": hit_last_str, "hit_pe_str": hit_pe_str, "ms_nf_str": ms_nf_str, "cs_nf_str": cs_nf_str, "me_nf_str": me_nf_str, "ce_nf_str": ce_nf_str, "shift_str": shift_str, "region1": region1, "region2": region2, "region3": region3}


    for sp_id in ir_ok0:
        w0_line = f'{sp_id}\tref_exon_ir\t{ir_ok0[sp_id]["tx_ir_str"]}\t{ir_ok0[sp_id]["tx_ok_str"]}\t{ir_ok0[sp_id]["region1"]}\t{ir_ok0[sp_id]["region2"]}\t{ir_ok0[sp_id]["region3"]}\n'
        if (w0_line in w0) == False:
            writer0.write(w0_line)
            w0.append(w0_line)

    for sp_id in ir_ok1:
        w1_line = f'{sp_id}\tref_exon_ir_poison_by_ir\t{ir_ok1[sp_id]["pe_len"]}\t{ir_ok1[sp_id]["avg_phylop"]}\t{ir_ok1[sp_id]["med_phylop"]}\t{ir_ok1[sp_id]["ptc_tx_str"]}\t{ir_ok1[sp_id]["shift_str"]}\t{ir_ok1[sp_id]["ms_nf_str"]}\t{ir_ok1[sp_id]["cs_nf_str"]}\t{ir_ok1[sp_id]["me_nf_str"]}\t{ir_ok1[sp_id]["ce_nf_str"]}\t{ir_ok1[sp_id]["d2tss_str"]}\t{ir_ok1[sp_id]["hit_tail_str"]}\t{ir_ok1[sp_id]["hit_last_str"]}\t{ir_ok1[sp_id]["hit_pe_str"]}\t{ir_ok1[sp_id]["region1"]}\t{ir_ok1[sp_id]["region2"]}\t{ir_ok1[sp_id]["region3"]}\n'
        if (w1_line in w1) == False:
            writer1.write(w1_line)
            w1.append(w1_line)

    return


def run(gtf, fasta, phylop, output_folder, version):
    o1 = open(f'{output_folder}/spliceIDs_ref.{version}.txt', "w")
    o1.write(f'spliceID\tType\tTxMain\tTxAlt\tCoord_R1\tCoord_R2\tCoord_R3\n')

    if fasta != None:
        o2 = open(f'{output_folder}/spliceIDs_ref.poison_exon.{version}.txt', "w")
        o2.write(f'spliceID\tType\tPElength\tPEavgPhyloP\tPEmedPhyloP\tTxCDtoNMD\tORFshift\tmRNA_start_NF\tcds_start_NF\tmRNA_end_NF\tcds_end_NF\tPTCtoTSS\tPTCinLast2Exons\tPTCinLastExon\tPTCinPE\tCoord_R1\tCoord_R2\tCoord_R3\n')
        if phylop != None:
            pybw = pyBigWig.open(f'{phylop}')
        else:
            pybw = None
    else:
        o2 = None
        pybw = None

    i = 0
    input_lines = []
    with open(gtf, "r") as f:
        for line in f:
            if line[0] == "#":
                continue
            ch, src, tp, start, end, null1, st, null2, info = line.rstrip().split("\t")
            if tp == "gene":
                if len(input_lines) != 0:
                    ch1, st1, gi1, gn1, intron_pairs, tx_dic, exon_dic = tandem_exons(input_lines)
                    judge_exon_skipping(fasta, ch1, st1, gi1, gn1, intron_pairs, tx_dic, o1, o2, pybw)
                    judge_exon_altLen(fasta, ch1, st1, gi1, gn1, intron_pairs, tx_dic, o1, o2, pybw)
                    judge_exon_ir(fasta, ch1, st1, gi1, gn1, exon_dic, tx_dic, o1, o2, pybw)
                    input_lines = []
                gene_eid = get_attribute(info, "gene_id")
                i += 1
                print(f'{i}: {gene_eid}')
            input_lines.append(line)

    #process last gene
    if len(input_lines) != 0:
        ch1, st1, gi1, gn1, intron_pairs, tx_dic, exon_dic = tandem_exons(input_lines)
        judge_exon_skipping(fasta, ch1, st1, gi1, gn1, intron_pairs, tx_dic, o1, o2, pybw)
        judge_exon_altLen(fasta, ch1, st1, gi1, gn1, intron_pairs, tx_dic, o1, o2, pybw)
        judge_exon_ir(fasta, ch1, st1, gi1, gn1, exon_dic, tx_dic, o1, o2, pybw)
    
    o1.close()

    if fasta != None:
        o2.close()

    return



parser = argparse.ArgumentParser(description = 'Build Splicing Events from GTF File')
parser.add_argument('gtf', help='An input GTF file to build splicing events.')
parser.add_argument('-o', '--output_dir', help = '(Optional) Output directory. Default: the same directory of the input GTF.')
parser.add_argument('-n', '--name', default = 'Pass', help = '(Optional) Version of the input GTF. E.g. GRCh38.100.')
parser.add_argument('--poison_exon', action = 'store_true', help = '(Optional) When used, splicing events for poison exons will be built in a separate file.')
parser.add_argument('--fasta', help = 'Only used when --poison_exon is active. Otherwise ignored. A genome fasta file that is paired with GTF.')
parser.add_argument('--phylop', help = 'Optional only used when --poison_exon is active. Otherwise ignored. A bigwig file for PhyloP scores across species.')
args = parser.parse_args()


if args.poison_exon:
    if args.fasta == None:
        sys.exit("Error: --fasta must be specified when --poison_exon is active.")

if args.output_dir == None:
    if "/" in args.gtf:
        out_dir = "/".join(args.output_dir.split("/")[:-1])
    else:
        out_dir = "."
else:
    out_dir = args.output_dir


run(args.gtf, args.fasta, args.phylop, out_dir, args.name)
