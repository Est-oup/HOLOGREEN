//#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Project dir setup
project_dir = "${workflow.projectDir}"


// Input folder setup
query_folder = "${project_dir}/data/query_contigs"
ref_folder = "${project_dir}/data/seq_references"
bact_DNA_pol = "${project_dir}/data/DNA_pol_bact_DB"

// Output folder setup
   	  //GENERAL
db_folder = "${project_dir}/out/databases"
db_folder_mmseqs_query = "${project_dir}/out/databases/mmseqs/query"
db_folder_mmseqs_ref = "${project_dir}/out/databases/mmseqs/ref"
log_folder = "${project_dir}/out/log"


   	  //ETAPE 1: ORF SEARCH
	  
    //miniprotsearch
//Prediction ORF
orf_folder = "${project_dir}/out/1_Prediction_orf"
orf_folder_log = "${orf_folder}/log"

    //Prodigalsearch
prodigal_folder = "${orf_folder}/prodigal"
orfpred_prodigal_folder = "${prodigal_folder}/orf"
orfpred_prodigal_folder_temp = "${prodigal_folder}/temp"


//Miniprotsearch
//initialisation
orf_miniprot_folder = "${orf_folder}/miniprot"
orf_init_miniprot_fasta = "${orf_miniprot_folder}/orf_init_fasta"


//1st Alignement

db_folder_miniprot_query = "${orf_miniprot_folder}/databases/query"
miniprot_alignment_folder = "${orf_miniprot_folder}/alignment/"
miniprot_alignment_temp_folder = "${miniprot_alignment_folder}/temp"
miniprot_info_folder = "${orf_miniprot_folder}/info"

// orf finded
orfpred_miniprot = "${orf_miniprot_folder}/orfpredminiprot_brut"

// extraction
orfpred_extracted_miniprot = "${orf_miniprot_folder}/extracted_results"
orfpred_extracted_miniprot_detail = "${orfpred_extracted_miniprot}/orf_by_prot"
orfpred_extracted_miniprot_file = "${orfpred_extracted_miniprot}/miniprot_orf.fasta"
orfpred_extracted_miniprot_log = "${orfpred_extracted_miniprot}/miniprot_orf_extraction.log"


    //general output of ORF reasearch
orf_pred_file = "${orf_folder}/unclustered_ORF"
orf_pred_file_derep  = "${orf_folder}/clustered_ORF"


   	  //ETAPE 2: ALIGNEMENT PROTEIQUE GENERAL

    //Alignment folder
alignement_general = "${project_dir}/out/2_General_alignment"
ref_combined ="${alignement_general}/ref_combined"
alignment_folder = "${alignement_general}/alignement"
alignment_temp_folder = "${alignement_general}/alignement_temp"
alignment_filt = "${alignement_general}/alignement_filt"
alignment_filt_marker = "${alignement_general}/alignement_filt_marker"
contigs_duplicated = "${alignement_general}/contigs_duplicated"
orf_ids = "${alignement_general}/orf_ids"



         //ETAPE 3: BACT DNA POL remove
DNA_pol =  "${project_dir}/out/3_DNA_pol_thrim"
DNA_pol_fasta = "${DNA_pol}/fasta/"
DNA_pol_fasta_db = "${DNA_pol}/db/"
DNA_pol_files = "${DNA_pol}/DNA_pol_ORF/"
DNA_pol_align = "${DNA_pol}/alignement/aln"
DNA_pol_align_temp = "${DNA_pol}/alignement/temp"
DNA_pol_align_alnfilt = "${DNA_pol}/alignement/aln_filt"
DNA_pol_align_alnfilt_ID = "${DNA_pol}/alignement/aln_filt_ID"
DNA_pol_align_final_aln = "${DNA_pol}/alignement/final_aln"
DNA_pol_align_final_aln_temp = "${DNA_pol}/alignement/final_aln_temp"





        //ETAPE 5: Final treatment of folder
alignment_final = "${project_dir}/out/5_out_final/alignments"
orf_codants_final = "${project_dir}/out/out_final/"

workflow {
	 //ETAPE 1: Recherche d'ORF
	 
    // Etape 1A: Prodigeal search
    log_1A = _1A_prodigalsearch()

    // Etape 1A: Initialisation recherche d'ORF avec alignement sur prot de séquences
    log_1B = _1B_miniprotinit(log_1A)

    // Etape 1A: Initialisation recherche d'ORF avec alignement sur prot de séquences
    log_1C = _1C_miniprotinit(log_1B)

    // Etape 1B: Recherche d'ORF potentiels using miniprot 
    log_1D = _1D_miniprotsearch(log_1C)

    // Etape 1E : extrcation des résultats
    log_1E = _1E_extract_miniprot_results(log_1D)

    // Etape 1D: Combined ORF
    log_1F = _1F_combinedorf(log_1E)

   	 //ETAPE 2: ALIGNEMENT PROTEIQUE GENERAL
    log_2 = _2_general_alignment(log_1F)

	 //ETAPE 3: BACT DNA POL remove	 
    log_3 = _3_DNApolthrim(log_2)

}



			// ETAPE 1 : PREDICTION ORF 

//Etape 1A:
process _1A_prodigalsearch {
    conda 'bioconda::prodigal conda-forge::biopython'

    output:
    path 'log_1A'

    script:
    """
    mkdir -p "${orfpred_prodigal_folder}"
    mkdir -p "${orfpred_prodigal_folder_temp}"
    mkdir -p "${log_folder}"
    python "${project_dir}/scripts/STEP_1/STEP1A_prodigal.py" ${query_folder} ${orfpred_prodigal_folder} ${orfpred_prodigal_folder_temp} ${prodigal_folder}/log_prodigal.txt

    #log
    echo "log_1A: processus 1A recherche d'ORF par prodigal : terminé" > "${log_folder}/log_STEP1"
    echo "log_1A" > log_1A
    """
}


process _1B_miniprotinit {
    conda 'bioconda::mmseqs2 conda-forge::biopython'

    input:
    path log_1A

    output:
    path 'log_1B'

    script:
    """ 
    #Extraction des contigs qui ont des ORF
    mkdir -p "${orfpred_prodigal_folder}"
    mkdir -p "${orf_folder_log}"
    mkdir -p ${orf_init_miniprot_fasta}
    python "${project_dir}/scripts/STEP_1/STEP1B_init_orf_fasta.py" ${orfpred_prodigal_folder} ${query_folder} ${orf_init_miniprot_fasta}/contigs_with_orf.fasta ${orf_folder_log}/log_1B
    #log
    echo "log_1B: extraction des sequences des ORF terminé" >> "${log_folder}/log_STEP1"

    #Création des db pour blast 
    mkdir -p "${db_folder_mmseqs_ref}"
    python "${project_dir}/scripts/STEP_1/STEP1B_create_db.py" ${ref_folder} ${db_folder_mmseqs_ref}
    #log
    echo "log_1B: création db references terminé" >> "${log_folder}/log_STEP1"

    #Blast des contigs avec les prot de ref pour voir les contigs qui ont un potentiel ORF
    mkdir -p "${miniprot_alignment_temp_folder}"
    python "${project_dir}/scripts/STEP_1/STEP1B_blast_mmseqs.py" ${orf_init_miniprot_fasta}/contigs_with_orf.fasta ${db_folder_mmseqs_ref} ${miniprot_alignment_folder} ${miniprot_alignment_temp_folder}
    
    #log
    echo "log_1B: blast des contigs avec refdatabase terminé" >> "${log_folder}/log_STEP1"
    #
    #log
    echo "log_1B: processus terminé" >> "${log_folder}/log_STEP1"
    echo "log_1B" > log_1B
    """
}


//Etape 1C:
process _1C_miniprotinit {
    conda 'conda-forge::biopython'

    input:
    path log_1B

    output:
    path 'log_1C'

    script:
    """   
    # Traitement des alignements
    mkdir -p "${miniprot_info_folder}"
    python "${project_dir}/scripts/STEP_1/STEP1C_process_alignments.py" ${miniprot_alignment_folder} ${miniprot_info_folder}
    
    #log
    echo "log_1C: traitement des alignements ID terminé" >> "${log_folder}/log_STEP1"
    echo "log_1C: processus terminé" >> "${log_folder}/log_STEP1"
    echo "log_1C" > log_1C
    """
} 


//Etape 1D:
process _1D_miniprotsearch {
    conda 'conda-forge::biopython'

    input:
    path log_1C

    output:
    path 'log_1D'

    script:
    """ 
    conda install bioconda::miniprot

    #log 
    echo "log_1D: environnement conda terminé" >> "${log_folder}/log_STEP1"
    
    # Execution du script Python pour la prédiction des ORF
    mkdir -p "${orfpred_miniprot}"
    python "${project_dir}/scripts/STEP_1/STEP1D_orf_prediction_miniprot.py" ${miniprot_info_folder}  ${query_folder} ${ref_folder} ${orfpred_miniprot} ${orf_folder_log}
    
    # Log
    echo "log_1D: Prédiction des ORF avec miniprot terminé" >> "${log_folder}/log_STEP1"
    echo "log_1D" > log_1D
    """
}


// Extraction des résultats de Miniprot
process _1E_extract_miniprot_results {
    conda 'conda-forge::biopython'

    input:
    path log_1D

    output:
    path 'log_1E'

    script:
    """
    # extraction des résultats miniprot
    mkdir -p "${orfpred_extracted_miniprot}"
    mkdir -p "${orfpred_extracted_miniprot_detail}"
    python "${project_dir}/scripts/STEP_1/STEP1E_extract_miniprot_results.py" ${orfpred_miniprot} ${orfpred_extracted_miniprot_detail} ${orfpred_extracted_miniprot_file} "${orf_folder_log}/log_1E"

    #log
    echo "log_1E: Extraction des résultats de Miniprot terminée" >> "${log_folder}/log_STEP1"
    echo "log_1E" > log_1E
    """
}




//Etape 1F:
process _1F_combinedorf {
    conda 'bioconda::cd-hit conda-forge::biopython'

    input:
    path log_1E

    output:
    path 'log_1F'

    script:
    """
    mkdir -p ${orf_pred_file}
    mkdir -p ${orf_pred_file_derep}
    mkdir -p ${orf_pred_file_derep}/logs

    python "${project_dir}/scripts/STEP_1/STEP1F_combined_orf.py" ${orfpred_prodigal_folder} ${orfpred_extracted_miniprot_detail} "${orf_pred_file}/ORF_pred" "${orf_pred_file_derep}/ORF_pred_derep" "${orf_pred_file_derep}/logs/ORF_pred_derep.log"
    
    #log
    echo "log_1F: processus 1F merge des deux prédictions d'ORF : terminé" >> "${log_folder}/log_STEP1"
    echo "log_1F" > log_1F
    """
}

		// ETAPE 2: ALIGNEMENT PROTEIQUE GENERAL
				
//ETAPE 2:
process _2_general_alignment {
    conda 'bioconda::mmseqs2 conda-forge::biopython'

    input:
    path log_1F

    output:
    path 'log_2'

    script:
    """
    #Etape 2A: Initialisation de l'alignement MMseqs2

    #Creation de la db orf 
    python "${project_dir}/scripts/STEP_2/STEP2A_create_db_mmseqs.py" "${orf_pred_file_derep}/ORF_pred_derep.final" ${db_folder_mmseqs_query}

    #Combinaisons des fichiers de reference
    python "${project_dir}/scripts/STEP_2/STEP2A_combine_ref.py" ${ref_folder} ${ref_combined}

    #log
    echo "log_2A: creation db query mmseqs : terminé" > "${log_folder}/log_STEP2"   


    #Etape 2B: Effectuer la recherche MMseqs2
    python "${project_dir}/scripts/STEP_2/STEP2B_blast_mmseqs.py" ${db_folder_mmseqs_query} ${ref_combined} ${alignment_folder} ${alignment_temp_folder}

    #log
    echo "log_2B: alignement general mmseqs : terminé" > "${log_folder}/log_STEP2"    

    #Etape 2C: Extraction informations alignement 
    python "${project_dir}/scripts/STEP_2/STEP2C_process_alignment.py" ${alignment_folder} ${alignment_filt} ${alignment_filt_marker} ${contigs_duplicated} ${orf_ids}

    #log
    echo "log_2C: extraction des informations : terminé" >> "${log_folder}/log_STEP2"

    #Etape 2D: Extraire les séquences FASTA à partir des identifiants
    python "${project_dir}/scripts/STEP_2/STEP2D_extractfastafromID.py" "${orf_pred_file_derep}/ORF_pred_derep.final" ${orf_ids}  ${orf_codants_final}

    #log
    echo "log_2D: extraction des seqs fasta : terminé" >> "${log_folder}/log_STEP2"
    echo "log_2: processus : terminé" >> "${log_folder}/log_STEP2"
    echo "log_2" > log_2 
    """
}





			// ETAPE 3 : BACT DNA POL remove 

//Etape 3:
process _3_DNApolthrim {
    conda 'bioconda::mmseqs2 conda-forge::pandas conda-forge::biopython'


    input:
    path log_2

    output:
    path 'log_3'
    
    script:
    """
    # Etape 3A: Initialisation recherche d'ORF avec alignement sur prot de séquences
    #Renomer dans un nouveau dossier par virus et par bacteria
    #Création des db references de virus et bactéries
    mkdir -p "${DNA_pol_fasta}"
    python "${project_dir}/scripts/STEP_3/STEP3A_DNA_pol_init.py" ${bact_DNA_pol} "${ref_folder}/dnapol.fasta" ${DNA_pol_fasta}
    #echo "log_3A: STEP3_DNA_pol_init : terminé" > "${log_folder}/log_STEP3"    
    
    #Etape 3B:Construction DB par mmseqs
    python "${project_dir}/scripts/STEP_3/STEP3B_create_db_mmseqs.py" ${DNA_pol_fasta} ${DNA_pol_fasta_db}
    echo "log_3B: Construction DB par mmseqs : terminé" >> "${log_folder}/log_STEP3"
 
    #Etape 3C: Alignement DB avec les ORF pred  
    python "${project_dir}/scripts/STEP_3/STEP3C_filter_file.py" ${orf_codants_final} ${DNA_pol_files}
    echo "log_3C: Alignement DB avec les ORF pred : terminé" >> "${log_folder}/log_STEP3"    

    #Etape 3D : blast 
    python "${project_dir}/scripts/STEP_3/STEP3D_blast_mmseqs.py" ${DNA_pol_fasta_db} ${DNA_pol_files} ${DNA_pol_align} ${DNA_pol_align_temp}
    echo "log_3D: blast : terminé" >> "${log_folder}/log_STEP3"

    #Etape 3E: receuil des infos
    python "${project_dir}/scripts/STEP_3/STEP3E_DNA_pol_filt.py" ${DNA_pol_align} ${DNA_pol_align_alnfilt} ${DNA_pol_align_alnfilt_ID}
    echo "log_3E: receuil des infos : terminé" >> "${log_folder}/log_STEP3"

    #Etape 3F: recuperation des fastas
    python "${project_dir}/scripts/STEP_3/STEP3F_extractfastafromID.py" "${orf_codants_final}/dnapol.fasta" ${DNA_pol_align_alnfilt_ID} ${orf_codants_final}
    echo "log_3F: recuperation des fastas : terminé" >> "${log_folder}/log_STEP3"

    #log
    echo "log_3" > log_3 
    """
}

