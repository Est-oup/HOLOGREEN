#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Project dir setup
project_dir = "${workflow.projectDir}"
scripts = "${project_dir}/scripts"
out = "${project_dir}/out"

// Input folder setup
orf_codants =  "${project_dir}/data/orf_codants"
ref_folder = "${project_dir}/data/seq_references"
outgroup_folder = "${project_dir}/data/outgroup"
clusters = "${project_dir}/data/bin/clusters.tsv"

// Output folder

//STEP 1:
step_1 = "${out}/STEP_1"
alignments_init = "${step_1}/alignment"
alignments_init_temp = "${step_1}/alignment_temp"
db_folder_mmseqs_ref = "${step_1}/database"


//STEP 2:
step_2 = "${out}/STEP_2"
ref_db = "${step_2}/ref_db"
bounds_alignment = "${step_2}/bounds_alignment"
bounds_alignment_temp = "${step_2}/bounds_alignment_temp"
bounds_results = "${step_2}/bounds_results"
bounds_summary = "${step_2}/bounds_summary"
filt_results = "${step_2}/filt_results"
score_aln_filt = "${step_2}/score_aln_filt"
orf_codants1 = "${step_2}/orf_codants_score"
score_filt = "${step_2}/filt_results/score_filt"
log_analysis = "${out}/log_analysis"
bounds_analysis = "${log_analysis}/bounds_performance"

//ETAPE 3: ORF filtration by size and clustered
step_3 = "${out}/STEP_3"
minimum_size = "${step_3}/minimum_size.txt"
orf_codants2 = "${step_3}/orf_codant2" //orf filtered by size

//ETAPE 4: deduplication and keep best alignment
step_4 = "${out}/STEP_4"
orf_codants3 = "${step_4}/orf_codants3"
orf_codants4 = "${step_4}/orf_codants4"
best_alignment = "${step_4}/best_alignment"

//ETAPE 5: Performence analysis
general_analysis = "${log_analysis}/general_performance"

//ETAPE 6: Alignment
step_6 = "${out}/STEP_6"
init_seqs = "${step_6}/init_seqs"
alignment = "${step_6}/alignment"
alignment_filt = "${step_6}/alignment_filt"

//ETAPE 7: Fitration gap alignment et log
step_7 = "${out}/STEP_7"
tree_data = "${step_7}/tree_data"
tree_data_filt = "${step_7}/tree_data_filt"
tree_pdf = "${step_7}/tree_pdf"

// Log folder setup
log_folder = "${project_dir}/logs"



workflow {
    //ETAPE 1 : Alignement général 
    log_1 = _1_general_alignment()

    //ETAPE 2 : Filtration of ORF alignments based on evalue et reference score
    log_2 = _2_filter_eval_refbounds(log_1)

    // ETAPE 3 : Filtration by size and clustering
    log_3 = _3_size_filter(log_2)

    // ETAPE 4 : Remove duplicated ORF within contigs
    log_4 = _4_dup_remove(log_3)
    
    // ETAPE 4 : Analysis of filtration 
    log_5 = _5_analysis(log_4)

    //ETAPE 6 :
    log_6 =  _6_alignment(log_5)
    
    //ETAPE 7 : 
    log_7 = _7_tree(log_6)
}


// Etape 1 : Filter by evalue and ref bounds
process _1_general_alignment {
    conda 'bioconda::mmseqs2 conda-forge::biopython'

    output:
    path 'log_1'

    script:
    """
    #Création des db pour blast 
    mkdir -p "${db_folder_mmseqs_ref}"
    python "${project_dir}/scripts/STEP1A_create_db_mmseqs.py" ${ref_folder} ${db_folder_mmseqs_ref}
    #log
    echo "log_1: création db references terminé" >> "${log_folder}/log_STEP1"

    #Blast des contigs avec les prot de ref pour voir les contigs qui ont un potentiel ORF
    mkdir -p "${alignments_init}"
    mkdir -p "${alignments_init_temp}"
    python "${project_dir}/scripts/STEP1B_blast_mmseqs.py" ${orf_codants} ${db_folder_mmseqs_ref} ${alignments_init} ${alignments_init_temp}
    
    #log
    echo "log_1: blast des contigs avec refdatabase terminé" >> "${log_folder}/log_STEP1"
    #
    #log
    echo "log_1: processus terminé" >> "${log_folder}/log_STEP1"
    echo "log_1" > log_1
    """
}


// Etape 2 : Filter by evalue and ref bounds
process _2_filter_eval_refbounds {
    conda 'bioconda::mmseqs2 conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn conda-forge::biopython conda-forge::pyqt conda-forge::qt'

    input:
    path log_1

    output:
    path 'log_2'

    script:
    """
    #Setup out folder    
    mkdir -p ${step_2}
    mkdir -p ${score_filt}

    #Filtration by evalue score in alignments
    python ${scripts}/STEP2A_filt_eval.py ${alignments_init} ${filt_results} "1e-0,1e-5,1e-10,1e-20,1e-50,1e-75,1e-100,1e-150,1e-200"
        #Log
	echo "log_2: filtration of ORF by evalue : completed" >> "${log_folder}/log_STEP2"

    #Determine the score bounds of alignment for each marker gene
    python  ${scripts}/STEP2B_define_bounds.py ${ref_folder} ${ref_db} ${bounds_alignment} ${bounds_results} ${bounds_alignment_temp}
        #Log
        echo "log_2: determined bounds of ref : completed" >> "${log_folder}/log_STEP2"

    #Filter by the right bounds threashold for ORF
    python ${scripts}/STEP2C_filter_by_ref_score.py ${bounds_results} ${alignments_init} ${score_aln_filt} ${score_filt} ${orf_codants} ${orf_codants1} ${bounds_summary}
        #Log
	echo "log_2: filtration of ORF by socre bounds : completed" >> "${log_folder}/log_STEP2"

    #Filtration step analysis 
    mkdir -p ${bounds_analysis}
    python ${scripts}/STEP2D_analyze_eval_bounds.py ${filt_results} ${bounds_analysis} ${orf_codants} ${orf_codants1}

    #Log
    echo "log_2: Analysis of different filtration step : completed" >> "${log_folder}/log_STEP2"
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
    mkdir -p ${step_3}

    #Detection of the right size by marker
    python  ${scripts}/STEP3A_size_references.py ${ref_folder} ${minimum_size}
        #Log
        echo "log_3: size of referece detection : completed" > "${log_folder}/log_STEP3"

    #Filter the orf by the size
    python ${scripts}/STEP3B_filter_size.py ${orf_codants1} ${minimum_size} ${orf_codants2}


    #Log
    echo "log_3: filtration of ORFs size : completed" >> "${log_folder}/log_STEP3"
    echo "log_3" > log_3
    """
}

// Etape 4 : Remove duplication
process _4_dup_remove {
    conda 'conda-forge::biopython  conda-forge::pandas'

    input:
    path log_3

    output:
    path 'log_4'

    script:
    """
    #Setup out folder    
    mkdir -p ${step_4}

    #Filtration of the duplicated ORF by contigs and keep ORF with the best alignment
    python  ${scripts}/STEP4_filter_orf_dup.py ${alignments_init} ${orf_codants2} ${clusters} ${best_alignment} ${orf_codants3} ${orf_codants4}

    #Log
    echo "log_4: Deduplication of ORF : completed" > "${log_folder}/log_STEP4"
    echo "log_4" > log_4
    """
}

// Etape 5 : Analysis of ORF filtration
process _5_analysis {
    conda 'bioconda::mmseqs2 conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn conda-forge::biopython conda-forge::pyqt conda-forge::qt'
    
    input:
    path log_4

    output:
    path 'log_5'

    script:
    """
    mkdir -p ${general_analysis}

    #Detection of the right size by marker
    python ${scripts}/STEP5_analysis_filt.py ${orf_codants} ${orf_codants1} ${orf_codants2} ${orf_codants3} ${orf_codants4} ${general_analysis}

    #Log
    echo "log_5: analysis of filtration: completed" > "${log_folder}/log_STEP5"
    echo "log_5" > log_5
    """
}


// Etape 6 : Sequence alignment
process _6_alignment {
    conda 'conda-forge::numpy conda-forge::biopython bioconda::mafft'

    input:
    path log_5

    output:
    path 'log_6'

    script:
    """
    mkdir -p ${step_6}

    #Sequence file preparation (orf + references sequence concatenantion)
    python  ${scripts}/STEP6A_init_fasta.py ${orf_codants4} ${ref_folder} ${outgroup_folder} ${init_seqs}
    	        #Log
		echo "log_6A: File preparation : completed" > "${log_folder}/log_STEP6"

    #1st alignment 
    python  ${scripts}/STEP6B_align_sequences.py ${init_seqs} ${alignment}
    	        #Log
		echo "log_6B: alignment : completed" >> "${log_folder}/log_STEP6"
		
    #Alignment filtration  : position with 90% of gap is removed
    python  ${scripts}/STEP6C_filter_alignments.py ${alignment} ${alignment_filt} 90 ${log_analysis}/alignment_stats.csv
    	        #Log
		echo "log_6C: alignment filtration : completed" >> "${log_folder}/log_STEP6"

    #Log
    echo "log_6: Alignment step : completed" >> "${log_folder}/log_STEP6"
    echo "log_6" > log_6
    """
}

// Etape 7 : Tree calculation
process _7_tree {
    conda 'conda-forge::ete3 conda-forge::biopython bioconda::iqtree'
    
    input:
    path log_6

    output:
    path 'log_7'

    script:
    """
    mkdir -p ${step_7}


    #calculate tree
    python  ${scripts}/STEP7A_calculate_tree.py ${alignment_filt} ${tree_data}
    	        #Log
		echo "log_7A: calculate tree : completed" > "${log_folder}/log_STEP7"

    #Tree labels filtrations
    python  ${scripts}/STEP7B_newick_filt.py ${tree_data} ${ref_folder} ${tree_data_filt}
    	        #Log
		echo "log_7B: labels modifications : completed" >> "${log_folder}/log_STEP7"

    #Generate pdf 
    python  ${scripts}/STEP7C_generate_pdf.py ${tree_data_filt} ${tree_pdf} ${ref_folder}
    	        #Log
		echo "log_7C: pdf generation : completed" >> "${log_folder}/log_STEP7"

    #Log
    echo "log_7: Deduplication of ORF : completed" >> "${log_folder}/log_STEP7"
    echo "log_7" > log_7
    """
}