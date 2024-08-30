#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.time.LocalDateTime
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 

process RunNeuralNetwork{
    label "modidec"
    publishDir (path: "${params.out_dir}/", mode: "copy")
    stageInMode "copy"
    input:
        path reference_path
        path pod5_path
        path bam_path
        path model_path
        path level_table_file
    output:
        path("*.html")
    script:
    """
    mkdir -p pod5s
    mv *.pod5 pod5s/
    mv \$(basename ${params.model_path}) model
    python ${projectDir}/bin/analysis_neural_network.py -s ${params.start_index} -e ${params.end_index} -c ${params.chunk_size} -x ${params.max_seq_length} -r $reference_path -p ./pod5s -b $bam_path -m ./model -l $level_table_file 
    """
}

workflow {
    input_dirs = Channel.fromPath("${params.model_path}", type: 'dir')
    RunNeuralNetwork(file("${params.reference_path}"),file("${params.pod5_path}/*.pod5"),file("${params.bam_path}"),input_dirs,file("${params.level_table_file}"))
}


