#!/usr/local/bin/python3

import pandas as pd
from dataclasses import dataclass
import argparse
import math
#from os.path import exists

driver_type ='hom_del'
@dataclass
class SequenceRange:
    """Class for the start and end of a range."""
    name: str
    transcript: str
    start: int
    end: int
    chrom: str
    total_cn: int
    def overlaps(self, other: "SequenceRange") -> bool:
        if self.chrom != other.chrom:
            return False
        return (other.start <= self.start <= other.end) or (other.start <= self.end <= other.end) or (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end)

my_parser = argparse.ArgumentParser(description='find the missing data in cnv files')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-ploidy',
                       type=float,
                       help='ploidy')
my_parser.add_argument('-gene_df',
                       type=str,
                       help='path to the df of genes and cooridnates of canonical transcript')
my_parser.add_argument('-somatic_cnv_vcf',
                       type=str,
                       help='path to the mtr input cnv of the sample')
my_parser.add_argument('-organ',
                       type=str,
                       help='organ')
####data input

args = my_parser.parse_args()
sample = args.sample
organ = args.organ
ploidy = args.ploidy
gene_df_path = args.gene_df
cnv_path = args.somatic_cnv_vcf
    
    
amps = list()
missing_gene_data_sample= list()
missing_data_genes_next_to_amps = list()
gene_df = pd.read_csv(gene_df_path) 
#gene_df=gene_df.dropna()

#set threshold for amplifications
if driver_type == 'amp':
    if ploidy <2.5:
        amp_threshold = 5
    elif ploidy >= 2.5:
        amp_threshold = 9
if driver_type == 'hom_del':
    amp_threshold = 0

    
#file_exists = exists(cnv_path):
try:    
    cnv = pd.read_csv(cnv_path, '\t')
    cnv['total_cn'] = cnv['major_cn'] + cnv['minor_cn']
    cnv['width'] = cnv['end'] - cnv['start']
    width = list(cnv['width'])
    total_cn = list(cnv['total_cn'])
    cnv['id'] = cnv['seqnames'].astype(str) + '_' + cnv['start'].astype(str) + '_' + cnv['end'].astype(str) + '_' + cnv['total_cn'].astype(str) +'_' + sample
    id_list = list(cnv['id'])
    for contig in range(len(total_cn)):
        if driver_type == 'amp':
            if total_cn[contig] >= amp_threshold: 
                amps.append(id_list[contig]) 
        if driver_type == 'hom_del':
            #if width[contig] < 21500000:
            #if width[contig] < 3000000 and width[contig] > 600000:
            #if width[contig] < 1500000:
            #if width[contig] < 1000000:
            #if width[contig] < 1500000 and width[contig] > 600000:
            #if width[contig] > 200000:
            if width[contig] < 1500000:
                if total_cn[contig] == amp_threshold: 
                    amps.append(id_list[contig])         
    #take the list of amps obtained in for loop above and convert to a table
    if len(amps) >0:
        amps_df = pd.DataFrame(amps)
        amps_df[[ 'chr', 'start', 'end','total_cn', 'sample']] = amps_df[0].str.split('_', 4, expand=True)
        amps_df.drop(columns=[0])
    else:
        amps_df = pd.DataFrame(columns=[0])
    ##for each contig (no matter if it is amplified) report whether any contig overlaps with gene of interest - this is to identify samples with no data for the gene of interest for i in range(len(gene_df.index)): ##need file with all coding gene name chromosome coordinates #for gene in gene_df:
    genes_in_amps = [[] for _ in range(len(amps_df.index))]
    for i in range(len(gene_df.index)):  
        gene  = SequenceRange(gene_df['gene_name'][i], gene_df['transcript_ID'][i], gene_df['start'][i], gene_df['end'][i], gene_df['chr'][i], 'total_cn_placeholder')
        ##find the genes with missing data
        contig_overlapping_gene = list()
        cnv_chr = cnv.loc[cnv['seqnames'].astype('str') == gene.chrom]
        cnv_chr.index = pd.RangeIndex(len(cnv_chr.index))
        contigs_after_gene = list()
        contigs_before_gene = list()
        for contig in range(len(list(cnv_chr['total_cn']))):
            contig_range= SequenceRange('place_holder', 'place_holder', int(math.floor(float(cnv_chr['start'][contig]))), int(math.floor(float(cnv_chr['end'][contig]))), str(cnv_chr['seqnames'][contig]), cnv_chr['total_cn'][contig])
            if contig_range.overlaps(gene):
                contig_overlapping_gene.append(id_list[contig])
            else:
                if contig_range.start > gene.end:
                    distance_from_gene = gene.end - contig_range.start
                    contig_id = contig_range.chrom + '_' +str(contig_range.start) + '_' +str(contig_range.end) + '_' + str(contig_range.total_cn) + '_' + str(distance_from_gene) + '_' +sample
                    contigs_after_gene.append(contig_id)
                elif contig_range.end < gene.start:
                    distance_from_gene = gene.start - contig_range.end
                    contig_id = contig_range.chrom + '_' +str(contig_range.start) + '_' +str(contig_range.end) + '_' + str(contig_range.total_cn)+ '_'+ str(distance_from_gene) +'_' +sample
                    contigs_before_gene.append(contig_id)
        if len(contig_overlapping_gene)== 0:
            missing_gene_data_sample.append(gene.name + '_' + gene.transcript + '_' +str(gene.start) + '_' +str(gene.end) +'_' +gene.chrom + '_' + sample)
            
            ##if there are contigs after the missing gene. Find the closest one and see if it is amplified
            total_cn_of_amp_neighbours = []
            if len(contigs_after_gene) >0:
                contigs_after_gene_df = pd.DataFrame(contigs_after_gene)
                contigs_after_gene_df[[ 'chr', 'start', 'end','total_cn', 'distance_from_gene', 'sample']] = contigs_after_gene_df[0].str.split('_', 5, expand=True)
                contigs_after_gene_df['total_cn'] = contigs_after_gene_df['total_cn'].astype('int') 
                contigs_after_gene_df['distance_from_gene'] = contigs_after_gene_df['distance_from_gene'].astype('int') 
                
                if driver_type == 'amp':
                    if contigs_after_gene_df['total_cn'][contigs_after_gene_df['distance_from_gene'].idxmax()] > amp_threshold:
                        total_cn_of_amp_neighbours.append(contigs_after_gene_df['total_cn'][contigs_after_gene_df['distance_from_gene'].idxmax()])
                if driver_type == 'hom_del':   
                    if contigs_after_gene_df['total_cn'][contigs_after_gene_df['distance_from_gene'].idxmax()] == amp_threshold:
                        total_cn_of_amp_neighbours.append(contigs_after_gene_df['total_cn'][contigs_after_gene_df['distance_from_gene'].idxmax()])
            
            if len(contigs_before_gene) >0:    
                contigs_before_gene_df = pd.DataFrame(contigs_before_gene)
                contigs_before_gene_df[[ 'chr', 'start', 'end','total_cn', 'distance_from_gene', 'sample']] = contigs_before_gene_df[0].str.split('_', 5, expand=True)
                contigs_before_gene_df['total_cn'] = contigs_before_gene_df['total_cn'].astype('int') 
                contigs_before_gene_df['distance_from_gene'] = contigs_before_gene_df['distance_from_gene'].astype('int') 

                if driver_type == 'amp':
                    if contigs_before_gene_df['total_cn'][contigs_before_gene_df['distance_from_gene'].idxmin()] > amp_threshold:
                        total_cn_of_amp_neighbours.append(contigs_before_gene_df['total_cn'][contigs_before_gene_df['distance_from_gene'].idxmin()])
                if driver_type == 'hom_del':
                     if contigs_before_gene_df['total_cn'][contigs_before_gene_df['distance_from_gene'].idxmin()] == amp_threshold:
                        total_cn_of_amp_neighbours.append(contigs_before_gene_df['total_cn'][contigs_before_gene_df['distance_from_gene'].idxmin()])

                if len(total_cn_of_amp_neighbours) > 0:
                        missing_data_genes_next_to_amps.append(gene.name + '_' + gene.transcript + '_' +str(gene.start) + '_' +str(gene.end) +'_' +gene.chrom + '_' +str(total_cn_of_amp_neighbours) + '_' +sample)

        ##report the genes in each amp
        if len(amps) > 0:    
            for amp in range(len(amps_df.index)):
                amp_range= SequenceRange('place_holder', 'place_holder', int(math.floor(float(amps_df['start'][amp]))), int(math.floor(float(amps_df['end'][amp]))), str(amps_df['chr'][amp]), amps_df['total_cn'][amp])
                if amp_range.overlaps(gene) and amp_range.chrom == gene.chrom:
                    genes_in_amps[amp].append(gene.name +'_' + gene.transcript + '_' +str(gene.start) + '_' +str(gene.end) + '_' +gene.chrom + '_' + sample)
                    #genes_in_amps[amp].append(gene.name)
    amps_df['genes_in_amps'] = genes_in_amps
        
    #print(sample)
    with open(sample +  '_' +organ + '_mtr_format_cnv_missing.txt', 'w') as f:
        f.write(sample+' complete')
except FileNotFoundError as e:
    with open(sample + '_' +organ + '_mtr_format_cnv_missing.txt', 'w') as f:
        f.write(sample+' no mtr format cnv file')
        
if len(missing_gene_data_sample) >0:
    missing_data_samples_gene_df = pd.DataFrame(missing_gene_data_sample)
    #missing_data_samples_gene_df = missing_data_samples_mdm2_df.rename(columns={0: 'missing_mdm2_samples'})
    missing_data_samples_gene_df[['gene', 'transcript_ID', 'start', 'end','chr', 'sample']] = missing_data_samples_gene_df[0].str.split('_', 5, expand=True)
    missing_data_samples_gene_df.drop(columns=[0])
else:
    missing_data_samples_gene_df = pd.DataFrame(columns=[0])
    
if len(missing_data_genes_next_to_amps) >0:
    missing_data_genes_next_to_amps_df = pd.DataFrame(missing_data_genes_next_to_amps)
    #missing_data_samples_gene_df = missing_data_samples_mdm2_df.rename(columns={0: 'missing_mdm2_samples'})
    #missing_data_genes_next_to_amps_df[['gene', 'transcript_ID', 'start', 'end','chr', 'sample']] = missing_data_genes_next_to_amps_df[0].str.split('_', 5, expand=True)
    missing_data_genes_next_to_amps_df[['gene', 'transcript_ID', 'start', 'end','chr', 'total_cn_contig_after_contig_before', 'sample']] = missing_data_genes_next_to_amps_df[0].str.split('_', 6, expand=True)
    missing_data_genes_next_to_amps_df.drop(columns=[0])
else:
    missing_data_genes_next_to_amps_df = pd.DataFrame(columns=[0])

#output table of genes with missing data                                                                                                 
missing_data_samples_gene_df.to_csv(sample + '_' +organ +'_genes_with_missing_data.csv')

#output table of genes with missing data next to amps                                                                                          
missing_data_genes_next_to_amps_df.to_csv(sample +'_' +organ + '_genes_with_missing_data_next_to_hom_dels.csv')

#output amps_df
amps_df.to_csv(sample + '_' +organ + '_hom_dels.csv')
