#! /usr/bin/env nextflow

nextflow.enable.dsl=2


/**
===============================
Filter with kaiju pipeline
===============================

This Pipeline uses taxonomic classification by Kaiju (Menzel, P. et al.) to filter input fastq.gz sequences.

### Homepage / git
git@github.com:ikmb/filter_with_kaiju.git

**/

// Pipeline version

params.version = workflow.manifest.version


// Help message
helpMessage = """
  =================================================================
   IKMB | Filter with kaiju | version ${params.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  
  nextflow run filter.kaiju --samples sample.sheet.csv --kaiju_db /path/to/kaiju_db --keep U --library pe
  
  Mandatory arguments:
  --samples         CSV (comma-separated file) containing the information about samples' fastq.gz reads to be processed.
                    It can have as many columns as wanted. But first three columns must have headers as following:
                    sample,fastq_1,fastq_2,...
                    If single end reads, fill fastq_2 column with NA. Do not leave it empty.
                    
  --kaiju_db        The path to the directory containing kaiju database of files *.fmi and nodes.dmp.
   
  --keep            Reads processed by Kaiju to keep. Accepted values are: "C" (classified) or "U" (unclassified).
  
  --library         Indicate if reads are single-ed (value: "se") or paired-end (value: "pe").
  
  Optonal arguments:
  --outdir          Local directory to which all output is written (default: results).
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}


log.info "=========================================" 
log.info "=========================================" 
log.info "FILTER WITH KAIJU   P I P E L I N E"
log.info "IKMB pipeline version v${workflow.manifest.version}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "=== Inputs =============================="
log.info "Samples file:        ${params.samples}"
log.info "Kaiju database:      ${params.kaiju_db}"
log.info "Reads to keep:       ${params.keep}"
log.info "Sequencing library:  ${params.library}"
log.info "Output directory:    ${params.outdir}"
log.info "=========================================="
log.info "Command Line:         $workflow.commandLine"
if (workflow.containerEngine) {
	log.info "Container Engine: ${workflow.containerEngine}"
}
log.info "=========================================" 


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}


// Sanity checks 

    // Check if kaiju database is there
if (params.kaiju_db) {
   kaiju_db_fmi = file(params.kaiju_db + "/kaiju_db_*.fmi")
   kaiju_db_nodes = file(params.kaiju_db + "/nodes.dmp")
	
    if (!kaiju_db_nodes.exists()) {
            exit 1, "Could not find your kaiju DB nodes.dmp file - please check the path"
        } else if (kaiju_db_fmi) {
            if (!kaiju_db_fmi[0].exists()){
                exit 1, "Could not find your kaiju DB .fmi file - please check the path"
        }
    } else {
        exit 1, "Could not find your kaiju DB .fmi file - please check the path"}
} else {
 exit 1, "No kaiju database was specified, aborting..."

 }


// Define processes


process quality_check {

    label 'fastp'

	publishDir "${params.outdir}/fastp", mode: 'copy'
	
	conda 'bioconda::fastp=0.12.4'

	scratch true 
	
	
    input:
    tuple val(sample), file(fastq_1), file(fastq_2)

        
    output:
     tuple val(sample), path(reads)
        
    script:
    out_html = sample + ".html"
    out_json = sample + ".json"
    reads = sample + "*.qc.R{1,2}.fq.gz"
        
        if (params.library == "pe"){
            r1 = sample + ".qc.R1.fq.gz"
            r2 = sample + ".qc.R2.fq.gz"
            
            """
            fastp \
            -i $fastq_1 \
            -I $fastq_2 \
            -o $r1 \
            -O $r2 \
            -h $out_html \
            -j $out_json \
            -w ${task.cpus}
            """
        } else if (params.library == "se"){
            r1 = sample + ".qc.R1.fq.gz"
            
            """
            fastp \
            -i $fastq_1 \
            -o $r1 \
            -h $out_html \
            -j $out_json \
            -w ${task.cpus}
            """
        }    
}


process classify {
    label 'kaiju'
    
    conda 'bioconda::kaiju=1.8.2 conda-forge::libstdcxx-ng=11.2.0'

	scratch true 
	
    input:
    tuple val(sample), path(reads)
    
    output: 
    tuple val(sample), file(kaiju_out_keep)
    
    
    script:
    kaiju_out = sample + ".kaiju.out"
    kaiju_out_keep = sample + ".kaiju.out.keep"
    
    
     if (params.library == "pe"){
         
        r1 = sample + "*.qc.R1.fq.gz"
        r2 = sample + "*.qc.R2.fq.gz"
                
        """
        kaiju \
        -t ${params.kaiju_db}/nodes.dmp \
        -f ${params.kaiju_db}/kaiju_db_*.fmi \
        -i $r1 \
        -j $r2 \
        -o $kaiju_out \
        -z ${task.cpus} \
        -v

        grep "^$params.keep" $kaiju_out > $kaiju_out_keep
        """ 

        }
    else if (params.library == "se"){
        
        r1 = sample + "*.qc.R1.fq.gz"
                
        """
        kaiju \
        -t ${params.kaiju_db}/nodes.dmp \
        -f ${params.kaiju_db}/kaiju_db_*.fmi \
        -i $r1 \
        -o $kaiju_out \
        -z ${task.cpus} \
        -v

        grep "^$params.keep" $kaiju_out > $kaiju_out_keep
        """ 
    }    
} 

process extract_sequences {
    label 'seqtk'

	publishDir "${params.outdir}/seqtk", mode: 'copy'

	conda 'bioconda::seqtk=1.3'
	
	scratch true 
		
    input:
    tuple val(sample), path(kaiju_out_keep), file(fastq_1), file(fastq_2)
    
    output:
    tuple val(sample), path(clean_reads)

    script:
    clean_reads = sample + "*.qc.keep.R{1,2}.fq.gz"
    
    if (params.library == "pe"){
        r1_keep_temp = sample + ".qc.keep.R1.fq"
        r2_keep_temp = sample + ".qc.keep.R2.fq"
        r1_keep = sample + ".qc.keep.R1.fq.gz"
        r2_keep = sample + ".qc.keep.R2.fq.gz"
        
        """ 
        # Prepare Kaiju output to filter
        cut -f2 $kaiju_out_keep > keep.list
        
        seqtk subseq $fastq_1 keep.list > \
        $r1_keep_temp
        
        seqtk subseq $fastq_2 keep.list > \
        $r2_keep_temp
        
        #gzip
        
        gzip $r1_keep_temp
        gzip $r2_keep_temp    
        """
    

    } else if (params.library == "se"){
      
        r1_keep_temp = sample + ".qc.keep.R1.fq"
        r1_keep = sample + ".qc.keep.R1.fq.gz"
        
        """
        # Prepare Kaiju output to filter
        cut -f2 $kaiju_out_keep > keep.list
        
        seqtk subseq $fastq_1 keep.list > \
        $r1_keep_temp
        
        #gzip
        
        gzip $r1_keep_temp
        """ 
    }  
}

process store_kaiju {

    publishDir "${params.outdir}/kaiju", mode: 'copy'
    
    label 'gzip'

	scratch true 

    input:
    tuple val(sample), path(kaiju_out_keep)
    
    output:
    path kaiju_gzip

    script:
    kaiju_gzip = sample + ".kaiju.out.keep.gz"
    """ 
    gzip -f $kaiju_out_keep
    """ 
}


// Start workflow
workflow {

    Channel
        .fromPath(params.samples)
        .splitCsv(header:true)
        .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
        .set { samples_ch}


    quality_check(samples_ch)
    classify(quality_check.out)
    extract_sequences(classify.out.join(samples_ch, by: [0]))
    store_kaiju(classify.out) 
}














