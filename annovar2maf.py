#!/usr/bin/env python

# File: annovar2maf.py
# Author: Anand Mayakonda [https://github.com/PoisonAlien]
# Created: July 19, 2023
# Description: This script converts annovar annotations and bcftools csq output to MAF

# MIT License
# Copyright (c) [2023] [Anand Mayakonda]

import argparse
import os.path
import re


def get_variant_type(ref, alt, vc):
    """
    Estimate Variany Type based on ref and alt alleles
    """
    variant_type_mappings = {
        "Frameshift_INDEL": "INS" if len(alt) > len(ref) else "DEL",
        "Inframe_INDEL": "INS" if len(alt) > len(ref) else "DEL",
        "Missense_Mutation": "SNP"
    }

    ref_alt = f"{ref}>{alt}"
    ref_alt_len = len(ref) + len(alt)
    ref_alt_diff = len(ref) - len(alt)

    if vc in variant_type_mappings:
        return variant_type_mappings[vc]
    elif ref_alt_diff < 0:
        return "INS"
    elif ref_alt_diff > 0:
        return "DEL"
    elif ref_alt in ["->A", "->C", "->T", "->G"]:
        return "INS"
    elif ref_alt in ["A>-", "C>-", "T>-", "G>-"]:
        return "DEL"
    elif ref_alt_len == 2:
        return "SNP"
    elif ref_alt_len == 4:
        return "DNP"
    elif ref_alt_len == 6:
        return "TNP"
    elif ref_alt_len > 6:
        return "ONP"
    else:
        return "NA"


def reformat_aachange(input_string):
    """
    Reformat bcftools csq aminio acid change info to standard HGVSp format
    Example: 
    55HC>55HG to HC55HG 
    253Y to Y253Y 
    215E>215* to E215*
    """
    
    spl = input_string.split(">")
    if len(spl) == 2:
        hgvsp = re.sub('[^A-Za-z]+', '', spl[0]) + re.sub('[^0-9]+', '', spl[0]) + re.sub('[0-9]+', '', spl[1])
    else:
        hgvsp = re.sub('[^A-Za-z]+', '', spl[0]) + re.sub('[^0-9]+', '', spl[0]) + re.sub('[^A-Za-z]+', '', spl[0])
        
    return hgvsp

def csq2maf(csq, tsb):
    
    """
    Takes bcftools csq formatted output and converts to maf
    bcftools example command:
    bcftools norm -f GRCh37.fa -m -both -Oz foo.vcf.gz | bcftools csq -c CSQ -f GRCh37.fa -g Homo_sapiens.GRCh37.82.gff3.gz -p a | bcftools +split-vep /dev/stdin -Oz -o foo.csq.vcf.gz -c - -s worst
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%gene\t%transcript\t%Consequence\t%amino_acid_change\t%dna_change\n' foo.csq.vcf.gz > foo.csq.tsv
    
    foo.csq.tsv would be the input
    """
    
    csq_to_vc_dict = {"synonymous": "Silent",
                      "missense": "Missense_Mutation",
                      "stop_lost": "Nonstop_Mutation",
                      "stop_gained": "Nonsense_Mutation",
                      "inframe_deletion": "In_Frame_Del",
                      "inframe_insertion": "In_Frame_Ins",
                      "frameshift": "INDEL",
                      "splice_acceptor": "Splice_Site",
                      "splice_donor": "Splice_Site",
                      "start_lost": "Translation_Start_Site",
                      "splice_region": "Splice_Region",
                      "stop_retained": "Silent",
                      "5_prime_utr": "5'UTR",
                      "3_prime_utr": "3'UTR",
                      "non_coding": "RNA",
                      "intron": "Intron",
                      "intergenic": "IGR",
                      "inframe_altering": "",
                      "coding_sequence": "Missense_Mutation",
                      "feature_elongation": "Targeted_Region",
                      "start_retained": "Silent",
                      ".": "IGR",
                      "NA": "NA"}

    with open(csq, 'r') as csq_file:
        for line in csq_file:
            line_spl = line.strip().split("\t")
            variant_classification = csq_to_vc_dict.get(line_spl[6], "NA")
            # Add Variant-type annotations based on the difference between ref and alt alleles
            variant_type = get_variant_type(
                line_spl[2], line_spl[3], variant_classification)

            # Reformat amino acid change to HGVSp convention
            if line_spl[7] == ".":
                aachange = "NA"
            else:
                aachange = reformat_aachange(line_spl[7])

            # Refgene in Unknonw if IGR or missing
            if line_spl[4] in ["NA", "NONE", "."]:
                refgene = "Unknown"
            else:
                refgene = line_spl[4]

            maf = [meta, refgene, line_spl[0], line_spl[1],
                   line_spl[1], variant_classification, variant_type, line_spl[2], line_spl[3], line_spl[5], line_spl[8], line_spl[7], aachange]
            print('\t'.join([str(x) for x in maf]))

def parse_record(line, col_idx, col_idx_non):
    """
    Convert an annovar line to MAF. Retain all other columns as is
    """
    # Annovar to MAF mappings (http://annovar.openbioinformatics.org/en/latest/user-guide/gene/)
    annovar_values = {
        'exonic': 'RNA',
        'splicing': 'Splice_Site',
        'UTR5': "5'UTR",
        'UTR3': "3'UTR",
        'intronic': 'Intron',
        'upstream': "5'Flank",
        'downstream': "3'Flank",
        'intergenic': 'IGR',
        'frameshift insertion': 'Frame_Shift_Ins',
        'frameshift deletion': 'Frame_Shift_Del',
        'frameshift block substitution': 'Frameshift_INDEL',
        'frameshift substitution': 'Frameshift_INDEL',
        'stopgain': 'Nonsense_Mutation',
        'stoploss': 'Nonstop_Mutation',
        'startloss': 'Translation_Start_Site',
        'startgain': 'Unknown',
        'nonframeshift insertion': 'In_Frame_Ins',
        'nonframeshift deletion': 'In_Frame_Del',
        'nonframeshift block substitution': 'Inframe_INDEL',
        'nonframeshift substitution': 'Inframe_INDEL',
        'nonsynonymous SNV': 'Missense_Mutation',
        'synonymous SNV': 'Silent',
        'unknown': 'Unknown',
        'ncRNA_exonic': 'RNA',
        'ncRNA_intronic': 'RNA',
        'ncRNA_UTR3': 'RNA',
        'ncRNA_UTR5': 'RNA',
        'ncRNA': 'RNA',
        'ncRNA_splicing': 'RNA',
        'NA': 'NA',
        '.': 'NA'
    }
    
    # Take the fisrt functional entry
    linespl = line.split("\t")
    func_refgene = linespl[col_idx.get("Func.refGene")]
    func_refgene = func_refgene.split(";")[0]
    
    #Take the fisrt functional entry
    exonicFunc_refgene = linespl[col_idx.get("ExonicFunc.refGene")]
    exonicFunc_refgene = exonicFunc_refgene.split(";")[0]
    
    # Refgene in Unknonw if IGR or missing
    refgene = linespl[col_idx.get("Gene.refGene")]
    refgene = refgene.split(";")[0]
    if refgene in ["NA", "NONE"]:
        refgene = "Unknown"
    
    #Use first transcript changes by default
    aa_change = linespl[col_idx.get("AAChange.refGene")]
    aa_change = aa_change.split(",")[0]
    aa_change = aa_change.split(":")
    
    # "Transcript_ID", "Exon_Number", "HGVSc", "HGVSp"
    if len(aa_change) == 5:
        aa_change = [aa_change[1], aa_change[2], aa_change[3], aa_change[4]]
    else:
        aa_change = ["NA", "NA", "NA", "NA"]
    
    if exonicFunc_refgene == "NA":
        variant_classification = annovar_values[func_refgene]
    else:
        variant_classification = annovar_values[exonicFunc_refgene]
        
    # Add Variant-type annotations based on the difference between ref and alt alleles
    variant_type = get_variant_type(linespl[col_idx.get("Ref")], linespl[col_idx.get("Alt")], variant_classification)
    
    res = [refgene, linespl[col_idx.get("Chr")], linespl[col_idx.get("Start")], linespl[col_idx.get("End")], 
           variant_classification, variant_type, linespl[col_idx.get("Ref")], linespl[col_idx.get("Alt")]] + aa_change
    
    for idx in list(col_idx_non.values()):
        if len(linespl)-1 >= idx:
            res.append(linespl[idx])
        else:
            res.append("NA")
    
    return '\t'.join([str(x) for x in res])


def Diff(li1, li2):
    """
    Tiny function to diff two lists. From internet, dont remeber where unfortunately :\
    """
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

def read_annovar_file(file_path, meta, protocol):
    with open(file_path, 'r') as annovar_file:
        first_line = next(annovar_file).strip()
        
        if protocol == "refGene":
            essential_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
                              'ExonicFunc.refGene', 'AAChange.refGene']
        else:
            essential_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene', 'GeneDetail.ensGene',
                              'ExonicFunc.ensGene', 'AAChange.ensGene']
        
        nonessential_cols = Diff(first_line.split("\t"), essential_cols)
        
    
        essential_cols_dict = {col: idx for idx, col in enumerate(first_line.split("\t")) if col in essential_cols}
        nonessential_cols_dict = {col: idx for idx, col in enumerate(first_line.split("\t")) if col in nonessential_cols}
        
        # In case ensGene is used as a protocol, rename the keys to refGene to harmonise the input
        if protocol == "ensGene":
            key_mapping = {'Chr': 'Chr', 'Start': 'Start', 'End': 'End', 'Func.ensGene': 'Func.refGene', 
                           'Gene.ensGene': 'Gene.refGene', 'GeneDetail.ensGene': 'GeneDetail.refGene',
                           'ExonicFunc.ensGene': 'ExonicFunc.refGene', 'AAChange.ensGene': 'AAChange.refGene'}
            for old_key, new_key in key_mapping.items():
                essential_cols_dict[new_key] = essential_cols_dict.pop(old_key)

        hdr = ["Tumor_Sample_Barcode", "NCBI_Build", "Center", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
               "Reference_Allele", "Tumor_Seq_Allele2","Transcript_ID", "Exon_Number", "HGVSc", "HGVSp"] + nonessential_cols

        yield '\t'.join(hdr)

        for line in annovar_file:
            maf_line = [meta, parse_record(line.strip(), essential_cols_dict, nonessential_cols_dict)]
            print('\t'.join([str(x) for x in maf_line]))
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert annovar and bcftools-csq annotations to MAF", prog="annovar2maf")
    parser.add_argument(
        'input', help="Annovar anotations file [Ex: myanno.hg19_multianno.txt] or a csq formatted file.")
    parser.add_argument("-t", "--tsb", help="Sample name. Default parses from the file name")
    parser.add_argument(
        "-b", "--build", help="Reference genome build [Default: hg38]", default="hg38")
    parser.add_argument(
        "-p", "--protocol", help="Protocol used to generate annovar annotations [Default: refGene]", default="refGene", choices=["refGene", "ensGene"])
    parser.add_argument(
        "-c", "--csq", help="Input file is a bcftools csq formatted output", action='store_true')

    args = parser.parse_args()

    if args.tsb is None:
        tsb = os.path.basename(args.input).split(".")[0]
    else:
        tsb = args.tsb
    
    meta = '\t'.join([str(x) for x in [tsb, args.build, "NA"]])
    
    if args.csq == True:
        hdr = ["Tumor_Sample_Barcode", "NCBI_Build", "Center", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
           "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "HGVSc", "amino_acid_change", "HGVSp"]
        print('\t'.join(hdr))
        
        csq2maf(args.input, tsb)
    else:
        for output_line in read_annovar_file(args.input, meta, args.protocol):
            print(output_line)
