#!/usr/local/bin/python3

import pandas as pd
from dataclasses import dataclass
import argparse
import math




##function to detect ovelaps between genomic segments e.g. overlap between a gene and an amplificaiont
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

#read in the arguments from nextflow
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
####data input

args = my_parser.parse_args()
sample = args.sample
ploidy = args.ploidy
gene_df_path = args.gene_df
cnv_path = args.somatic_cnv_vcf
    

amps = list()
homdels=list()


#set threshold for a contig to be considered amplificatied based on ploidy
if ploidy <2.5:
    amp_threshold = 5
elif ploidy >= 2.5:
    amp_threshold = 9

#set the threshold for considering a contiguous homozygously deleted to 0
homdel_threshold = 0

    
#read in the copy number file (expecting ascat style input, which are then modified below)
cnv = pd.read_table(cnv_path, sep='\t')

#read in the genes for analysis
gene_df = pd.read_csv(gene_df_path) 

##because checking if chromosome values match, prep the gene df by removing any chr, Chr strings before the chromosome number
gene_df['chr'] = gene_df['chr'].str.replace('chr','')
gene_df['chr'] = gene_df['chr'].str.replace('Chr','')


##prep the copy number file
cnv = cnv.rename(columns={'total.copy.number.inTumour':'total_cn','Chromosome':'seqnames','chromStart':'start','chromEnd':'end','minor.copy.number.inTumour':'minor_cn'})
cnv['major_cn'] = cnv['total_cn'] -cnv['minor_cn']
cnv['width'] = cnv['end'] - cnv['start']
width = list(cnv['width'])
total_cn = list(cnv['total_cn'])
cnv['id'] = cnv['seqnames'].astype(str) + '_' + cnv['start'].astype(str) + '_' + cnv['end'].astype(str) + '_' + cnv['total_cn'].astype(str) +'_' + sample
id_list = list(cnv['id'])

##identify the amplified or deleted contigs
for contig in range(len(total_cn)):
    ###amplification
        if total_cn[contig] >= amp_threshold: 
            amps.append(id_list[contig]) 
    ##homdel
        if total_cn[contig] == homdel_threshold : 
            homdels.append(id_list[contig])         
#take the list of amplification and homozygous deletions obtained in for loop above and convert to a table
if len(amps) >0:
    amps_df = pd.DataFrame(amps)
    amps_df[[ 'chr', 'start', 'end','total_cn', 'sample']] = amps_df[0].str.split('_', n=4, expand=True)
    amps_df.drop(columns=[0])
    amps_df['type']= 'AMP'

if len(homdels) >0:
    homdel_df = pd.DataFrame(homdels)
    homdel_df[[ 'chr', 'start', 'end','total_cn', 'sample']] = homdel_df[0].str.split('_', n=4, expand=True)
    homdel_df.drop(columns=[0])
    homdel_df['type']= 'HOMDEL'

if len(homdels) >0 and len(amps) >0:
    cn_events = pd.concat([amps_df, homdel_df], ignore_index=True)
elif  len(homdels) ==0 and len(amps) >0:
    cn_events = amps_df
elif len(homdels) >0 and len(amps) ==0:
    cn_events = homdel_df
else:
    cn_events = pd.DataFrame(columns=[0])



##annotated the CN events for overlap with a gene

genes_in_event = [[] for _ in range(len(cn_events.index))]
for i in range(len(gene_df.index)):  
    ##make the dataclass for the gene 
    gene  = SequenceRange(gene_df['gene_name'][i], gene_df['transcript_ID'][i], gene_df['start'][i], gene_df['end'][i], gene_df['chr'][i], 'total_cn_placeholder')
    ##report the genes which overlap each amplification or homozygous deletion event
    if len(cn_events.index) > 0:    
        for amp in range(len(cn_events.index)):
            amp_range= SequenceRange('place_holder', 'place_holder', int(cn_events['start'][amp]), int(cn_events['end'][amp]), str(cn_events['chr'][amp]), cn_events['total_cn'][amp])
            if amp_range.overlaps(gene) and amp_range.chrom == gene.chrom:
                genes_in_event[amp].append(gene.name +'_' + gene.transcript + '_' +str(gene.start) + '_' +str(gene.end) + '_' +gene.chrom + '_' + sample)
cn_events['genes_in_CN_event'] = genes_in_event

cn_events['genes_in_CN_event'] = cn_events['genes_in_CN_event'].astype('str') 

##make a friendlier output file
#d = {}
#for feature in set(cn_events['type']):
#    d[feature] = {}
#    tmp = cn_events[cn_events['type']==feature]
#    for gene in set(gene_df['gene_name']):
#        d[feature][gene] = {}
#        if tmp['genes_in_CN_event'].str.contains(gene+'_').any() ==True:
#            d[feature][gene] = feature
#        else:
#            d[feature][gene] = 'no_driver'

#output amps_df

#cn_events = pd.DataFrame.from_dict(d).T
cn_events['sample'] = sample
#cn_events['type'] = cn_events.index
cn_events.to_csv(sample + '_annotated_CN_events.csv',index=False)
