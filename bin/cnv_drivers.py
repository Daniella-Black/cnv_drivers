#!/usr/local/bin/python3
import pandas as pd
from dataclasses import dataclass
import argparse
#from os.path import exists

@dataclass
class SequenceRange:
    """Class for the start and end of a range."""
    name: str
    transcript: str
    start: int
    end: int
    chrom: str
    def overlaps(self, other: "SequenceRange") -> bool:
        if self.chrom != other.chrom:
            return False
        return (other.start <= self.start <= other.end) or (other.start <= self.end <= other.end) or (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end)

my_parser = argparse.ArgumentParser(description='find the missing data in cnv files')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-input_file',
                       type=str,
                       help='path to the germline vcf')


args = my_parser.parse_args()
sample = args.sample
input_file = args.input_file


apobec= SequenceRange('place_holder', 'place_holder', int(38952741), int(38992804), str('chr22'))

#read in the germline cnv file
cnv = pd.read_csv(input_file, comment='#', sep='\t', header=None)
##apobec gene is on chr22 so filter for that
cnv = cnv[cnv[0]=='chr22']
##filter for any value <DUP] and <DEL> indicating a polymorphism.
cnv = cnv[cnv[4]!='.']
##get the start and stop of the region
cnv[['DRAGEN','TYPE','region']] = cnv[2].str.split(':',expand=True)
cnv[['start','end']] = cnv['region'].str.split('-',expand=True)
cnv = cnv.reset_index(drop=True)
apobec_overlap = pd.DataFrame()
if len(cnv.index) > 0:    
    for amp in range(len(cnv.index)):
        cnv_range= SequenceRange('place_holder', 'place_holder', int(cnv['start'][amp]), int(cnv['end'][amp]), str(cnv[0][amp]))
        if cnv_range.overlaps(apobec) and cnv_range.chrom == apobec.chrom:
            apobec_overlap = pd.concat([apobec_overlap, cnv.iloc[[amp]] ])
            #genes_in_amps[amp].append(gene.name)
    apobec_overlap['sample']=sample 
    apobec_overlap = apobec_overlap.reset_index(drop=True)

#output amps_df
apobec_overlap.to_csv(sample + '_APOBEC3A_B_germline_polymorphism_overlap.csv',index=False)
