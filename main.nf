#! /usr/bin/env nextflow
nextflow.enable.dsl=1
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey,ploidy, organ,somatic_cnv_vcf, gene_df -> [tumour_sample_platekey,ploidy,organ, file(somatic_cnv_vcf),  file(gene_df)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    set val(tumour_sample_platekey),val(ploidy), val(organ), file(somatic_cnv_vcf),  file(gene_df) from ch_input

    output:
    file "*_mtr_format_cnv_missing.txt"
    file "*_genes_with_missing_data.csv"
    file "*_amplifications.csv"
 
    script:
    """
    cnv_drivers.py -sample '$tumour_sample_platekey' -somatic_cnv_vcf '$somatic_cnv_vcf' -ploidy '$ploidy' -gene_df '$gene_df' -organ '$organ'
    """ 
}
