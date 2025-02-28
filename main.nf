#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Launch dir'] = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir'] = workflow.projectDir
summary['User'] = workflow.userName
summary['Output dir'] = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/* --------------------
| Help Message |
--------------------- */

def helpMessage() {

log.info """
Usage:
The typical command for running the pipeline is as follows:

Mandatory:
--inputlist Input file


--max_time Maximum time (time unit)
(default: $params.max_time)

""".stripIndent()
}

project_dir = workflow.projectDir
run_date = new java.text.SimpleDateFormat("yyyy_MM_dd").format(new Date())



ch_input = Channel
.fromPath(params.inputlist)
.ifEmpty {
error "Cannot find input file: ${params.inputlist}"
}
.splitCsv(skip: 1) // Skip the header row
.map { row ->
tuple(row[0], file(row[1]), row[2], file(row[3])) // Create a tuple with three elements
}


process cnv_drivers {
tag { tumour_sample_platekey }

errorStrategy 'ignore'
input:
tuple val(tumour_sample_platekey), path(somatic_cnv_vcf), val(ploidy), path(gene_df)

output:
  file "*_annotated_CN_events.csv"
script:
"""
cnv_drivers.py -sample '$tumour_sample_platekey' -somatic_cnv_vcf '$somatic_cnv_vcf' -ploidy '$ploidy' -gene_df '$gene_df'
"""
}


workflow {

 drivers = cnv_drivers(ch_input).collectFile(name: 'annotated_CN_events.csv',cache:'lenient',keepHeader:true, storeDir: "${params.outdir}/call_drivers")

}
