#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{sample,input_file -> [sample, file(input_file)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    //publishDir "${params.outdir}/$sample", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    set val(sample), file(input_file) from ch_input

    output:
    file "*_APOBEC3A_B_germline_polymorphism_overlap.csv"
 
    script:
    """
    cnv_drivers.py -sample '$tumour_sample_platekey' -input_file '$input_file'
    """ 
}
