#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.time.LocalDateTime
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 


process load_kmer_tables {

    input:
    val(flowcell_type)
    output:
    path("*.txt"), emit:kmer_lvl_table
    script:

    //load kmer table from website
    """
    if [[ ${flowcell_type} == "RNA002" ]]; then
	wget https://raw.githubusercontent.com/nanoporetech/kmer_models/refs/heads/master/rna_r9.4_180mv_70bps/5mer_levels_v1.txt
	else
    wget https://raw.githubusercontent.com/nanoporetech/kmer_models/refs/heads/master/rna004/9mer_levels_v1.txt
    fi
    """
}

process RunNeuralNetwork{
    label "modidec"
    publishDir (path: "${params.out_dir}/", mode: "copy")
    stageInMode "symlink"
    input:
        path reference_path
        path pod5_path
        path bam_path
        path model_path
        path level_table_file
    output:
        path("*.html"), emit: html
    script:
    """
    mkdir -p pod5s
    mv *.pod5 pod5s/
    mv \$(basename ${params.model_path}) model
    python ${projectDir}/bin/analysis_neural_network.py -s ${params.start_index} -e ${params.end_index} -c ${params.chunk_size} -x ${params.max_seq_length} -r $reference_path -p ./pod5s -b $bam_path -m ./model -l $level_table_file -d ${params.mod_list} 
    """
}

workflow {
    input_dirs = Channel.fromPath("${params.model_path}", type: 'dir')
    kmer_table = load_kmer_tables(params.flowcell_type)
    RunNeuralNetwork(file("${params.reference_path}"),file("${params.pod5_path}/*.pod5"),file("${params.bam_path}"),input_dirs,kmer_table.kmer_lvl_table)
}


