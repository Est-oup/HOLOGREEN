#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Project dir setup
project_dir = "${workflow.projectDir}"
scripts = "${project_dir}/scripts"
out = "${project_dir}/out"

// Input folder setup
alignments_init = "${project_dir}/data/alignments"
orf_codants =  "${project_dir}/data/orf_codants"
ref_folder = "${project_dir}/data/seq_references"


// Output folder

//STEP 1:
step_1 = "${out}/STEP_1"
ref_db = "${step_1}/ref_db"
bounds_alignment = "${step_1}/bounds_alignment"
bounds_alignment_temp = "${step_1}/bounds_alignment_temp"
bounds_results = "${step_1}/bounds_results"
bounds_summary = "${step_1}/bounds_summary"
filt_results = "${step_1}/filt_results"
score_aln_filt = "${step_1}/score_aln_filt"
orf_codants1 = "${step_1}/orf_codants_score"
score_filt = "${step_1}/filt_results/score_filt"
log_analysis = "${out}/log_analysis"
bounds_analysis = "${log_analysis}/bounds_performance"


//ETAPE 2: deduplication and keep best alignment
step_2 = "${out}/STEP_2"
orf_codants_filt = "${step_2}/orf_codants_filt"
orf_codants2 = "${orf_codants_filt}/orf_codants2"
best_alignment = "${step_2}/best_alignment"

//ETAPE 3: ORF filtration by size and clustered
step_3 = "${out}/STEP_3"
minimum_size = "${step_3}/minimum_size.txt"
orf_codants3 = "${step_3}/orf_codants3" //orf filtered by size
orf_codants4 = "${step_3}/orf_codants4" // ORF clusterised final 
orf_codants5 = "${step_3}/orf_codants_clstr_info" // informations about orf clusters

//ETAPE 4: Performence analysis
general_analysis = "${log_analysis}/general_performance"

//ETAPE 5: Alignment
step_5 = "${out}/STEP_5"
init_seqs = "${step_5}/init_seqs"
alignment = "${step_5}/alignment"
alignment_filt = "${step_5}/alignment_filt"

//ETAPE 6: Fitration gap alignment et log
step_6 = "${out}/TSEP_6"
tree_data = "${step_6}/tree_data"
tree_data_filt = "${step_6}/tree_data_filt"
tree_pdf = "${step_6}/tree_pdf"

// Log folder setup
log_folder = "${project_dir}/logs"



workflow {
    //ETAPE 1 : Filtration of ORF alignments based on evalue et reference score
    log_1 = _1_filter_eval_refbounds()

    // ETAPE 2 : Remove duplicated ORF within contigs
    log_2 = _2_dup_remove(log_1)

    // ETAPE 3 : Filtration by size and clustering
    log_3 = _3_size_filter(log_2)

    // ETAPE 4 : Analysis of fitration 
    log_4 = _4_analysis(log_3)

    //ETAPE 5 :
    log_5 =  _5_alignment(log_4)
    
    //ETAPE 6 : 
    log_6 = _6_tree(log_5)
}


// Etape 1 : Filter by evalue and ref bounds
process _1_filter_eval_refbounds {
    conda 'bioconda::mmseqs2 conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn conda-forge::biopython'

    output:
    path 'log_1'

    script:
    """
    #Setup out folder    
    mkdir -p ${step_1}
    mkdir -p ${score_filt}

    #Filtration by evalue score in alignments
    python ${scripts}/STEP1A_filt_eval.py ${alignments_init} ${filt_results} "1e-0,1e-5,1e-10,1e-20,1e-50,1e-75,1e-100,1e-150,1e-200"
        #Log
	echo "log_1: filtration of ORF by evalue : completed" >> "${log_folder}/log_STEP1"

    #Determine the score bounds of alignment for each marker gene
    python  ${scripts}/STEP1B_define_bounds.py ${ref_folder} ${ref_db} ${bounds_alignment} ${bounds_results} ${bounds_alignment_temp}
        #Log
        echo "log_1: determined bounds of ref : completed" >> "${log_folder}/log_STEP1"

    #Filter by the right bounds threashold for ORF
    python ${scripts}/STEP1C_filter_by_ref_score.py ${bounds_results} ${alignments_init} ${score_aln_filt} ${score_filt} ${orf_codants} ${orf_codants1} ${bounds_summary}
        #Log
	echo "log_1: filtration of ORF by socre bounds : completed" >> "${log_folder}/log_STEP1"

    #Filtration step analysis 
    mkdir -p ${bounds_analysis}
    python ${scripts}/STEP1D_analyze_eval_bounds.py ${filt_results} ${bounds_analysis} ${orf_codants} ${orf_codants1}

    #Log
    echo "log_1: Analysis of different filtration step : completed" >> "${log_folder}/log_STEP1"
    echo "log_1" > log_1
    """
}

// Etape 2 : Remove duplication
process _2_dup_remove {
    conda 'conda-forge::biopython  conda-forge::pandas'

    input:
    path log_1

    output:
    path 'log_2'

    script:
    """
    #Setup out folder    
    mkdir -p ${step_2}
    mkdir -p ${orf_codants_filt}
    mkdir -p ${best_alignment}

    #Filtration of the duplicated ORF by contigs and keep ORF with the best alignment
    python  ${scripts}/STEP2_filter_orf_dup.py ${alignments_init} ${orf_codants} ${best_alignment} ${orf_codants2}

    #Log
    echo "log_2: Deduplication of ORF : completed" > "${log_folder}/log_STEP2"
    echo "log_2" > log_2
    """
}

// Etape 3 : filter ORF by size and clustering
process _3_size_filter {
    conda 'conda-forge::biopython  conda-forge::pandas'

    input:
    path log_2

    output:
    path 'log_3'

    script:
    """
    conda install bioconda::cd-hit

    mkdir -p ${step_3}

    #Detection of the right size by marker
    python  ${scripts}/STEP3A_size_references.py ${ref_folder} ${minimum_size}
        #Log
        echo "log_3: size of referece detection : completed" > "${log_folder}/log_STEP3"

    #Filtration of the orf by the size and clustering at 90%
    python ${scripts}/STEP3B_filter_size_and_cluster.py ${orf_codants2} ${minimum_size} ${orf_codants3} ${orf_codants4} ${orf_codants5} ${orf_codants5}/logs 0.90

    #Log
    echo "log_3: filtration of ORFs size and cluster: completed" >> "${log_folder}/log_STEP3"
    echo "log_3" > log_3
    """
}



// Etape 4 : Analysis of ORF filtration
process _4_analysis {
    conda 'bioconda::mmseqs2 conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn conda-forge::biopython'
    
    input:
    path log_3

    output:
    path 'log_4'

    script:
    """
    mkdir -p ${general_analysis}

    #Detection of the right size by marker
    python ${scripts}/STEP4_analysis_filt.py ${orf_codants} ${orf_codants1} ${orf_codants2} ${orf_codants3} ${orf_codants4} ${general_analysis}

    #Log
    echo "log_4: analysis of filtration: completed" > "${log_folder}/log_STEP4"
    echo "log_4" > log_4
    """
}


// Etape 5 : Sequence alignment
process _5_alignment {
    conda 'conda-forge::numpy conda-forge::biopython'

    input:
    path log_4

    output:
    path 'log_5'

    script:
    """
    mkdir -p ${step_5}

    #Sequence file preparation (orf + references sequence concatenantion)
    python  ${scripts}/STEP5A_init_fasta.py ${orf_codants4} ${ref_folder}  ${init_seqs}
    	        #Log
		echo "log_5A: File preparation : completed" > "${log_folder}/log_STEP5"

    conda install bioconda::mafft

    #1st alignment 
    python  ${scripts}/STEP5B_align_sequences.py ${init_seqs} ${alignment}
    	        #Log
		echo "log_5B: alignment : completed" >> "${log_folder}/log_STEP5"
		
    #Alignment filtration  : position with 90% of gap is removed
    python  ${scripts}/STEP5C_filter_alignments.py ${alignment} ${alignment_filt} 90 ${log_analysis}/alignment_stats.csv
    	        #Log
		echo "log_5C: alignment filtration : completed" >> "${log_folder}/log_STEP5"

    #Log
    echo "log_5: Alignment step : completed" >> "${log_folder}/log_STEP4"
    echo "log_5" > log_5
    """
}

// Etape 6 : Tree calculation
process _6_tree {
    conda 'conda-forge::ete3 conda-forge::biopython'
    
    input:
    path log_5

    output:
    path 'log_6'

    script:
    """
    conda install bioconda::fasttree

    mkdir -p ${step_6}


    #calculate tree
    python  ${scripts}/STEP6A_calculate_tree.py ${alignment_filt} ${tree_data}
    	        #Log
		echo "log_6A: calculate tree : completed" > "${log_folder}/log_STEP6"

    #Tree labels filtrations
    python  ${scripts}/STEP6B_newick_filt.py ${tree_data} ${ref_folder} ${tree_data_filt}
    	        #Log
		echo "log_6B: labels modifications : completed" >> "${log_folder}/log_STEP6"

    #Generate pdf 
    python  ${scripts}/STEP6C_generate_pdf.py ${tree_data_filt} ${tree_pdf} ${ref_folder}
    	        #Log
		echo "log_6C: pdf generation : completed" >> "${log_folder}/log_STEP6"

    #Log
    echo "log_6: Deduplication of ORF : completed" >> "${log_folder}/log_STEP6"
    echo "log_6" > log_6
    """
}