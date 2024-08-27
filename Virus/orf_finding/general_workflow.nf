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
orf_folder = "${project_dir}/out/pred_orf/"
orf_folder_log = "${orf_folder}/log"

    //Prodigalsearch
prodigal_folder = "${project_dir}/out/pred_orf/prodigal"
orfpred_prodigal_folder = "${prodigal_folder}/orf"
orfpred_prodigal_folder_temp = "${prodigal_folder}/temp"


//MIniprotsearch
//initialisation
orf_miniprot_folder = "${project_dir}/out/pred_orf/miniprot"
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
orf_pred_file = "${project_dir}/out/pred_orf/unclustered_ORF"
orf_pred_file_derep  = "${project_dir}/out/pred_orf/clustered_ORF"


   	  //ETAPE 2: ALIGNEMENT PROTEIQUE GENERAL

    //Alignment folder
alignment_folder = "${project_dir}/out/alignment"
alignment_temp_folder = "${alignment_folder}/temp"

    //Info folder
info_folder = "${alignment_folder}/info_alignment"
info_folder_IDS = "${info_folder}/raw"
info_folder_ID = "${info_folder}/dedup"
info_folder_countdup = "${info_folder}/countdup"
orf_codants_folder = "${info_folder}/orf_codants"


         //ETAPE 3: BACT DNA POL remove
DNA_pol =  "${project_dir}/out/DNA_pol_thrim"
DNA_pol_fasta = "${DNA_pol}/fasta/"
DNA_pol_fasta_db = "${DNA_pol}/db/"
DNA_pol_files = "${DNA_pol}/DNA_pol_ORF/"
DNA_pol_align = "${DNA_pol}/alignement/aln"
DNA_pol_align_temp = "${DNA_pol}/alignement/temp"
DNA_pol_align_alnfilt = "${DNA_pol}/alignement/aln_filt"
DNA_pol_align_alnfilt_ID = "${DNA_pol}/alignement/aln_filt_ID"
DNA_pol_align_seqfilt = "${DNA_pol}/alignement/seq_filt"
DNA_pol_align_final_aln = "${DNA_pol}/alignement/final_aln"
DNA_pol_align_final_aln_temp = "${DNA_pol}/alignement/final_aln_temp"



        //ETAPE 4: Polinton like thrim
polinto = "${project_dir}/out/polinto"
polinto_db = "${polinto}/db_ref"
polinto_fasta = "${polinto}/fasta_ref"
polinto_orf_init = "${polinto}/orf_init"
polinto_aln = "${polinto}/aln"
polinto_aln_temp = "${polinto}/aln_temp"
polinto_aln_filt = "${polinto}/aln_filt"
polinto_aln_filt_id = "${polinto}/aln_filt_id"
polinto_aln_final = "${polinto}/aln_final"
polinto_aln_final_temp = "${polinto}/aln_final_temp"
polinto_orf_final = "${polinto}/orf_final"


        //ETAPE 5: Final treatment of folder
alignment_final = "${project_dir}/out/out_final/alignments"
orf_codants_final = "${project_dir}/out/out_final/orf_codants"

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

    	 //ETAPE 3: BACT DNA POL remove	 
    log_4 = _4_polinto_thrim(log_3)

	 //ETAPE 5: Files results treatment	 
    log_5 = _5_final_treatment(log_4)


}



			// ETAPE 1 : PREDICTION ORF 

//Etape 1A:
process _1A_prodigalsearch {
    conda 'bioconda::prodigal'

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
    conda 'bioconda::mmseqs2'

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

    input:
    path log_1B

    output:
    path 'log_1C'

    script:
    """   
    # Traitement des alignements
    mkdir -p "${miniprot_info_folder}"
    python "${project_dir}/scripts/STEP_1/STEP1C_process_alignments.py" ${miniprot_alignment_folder} ${miniprot_info_folder}
    #
    #log
    echo "log_1C: traitement des alignements ID terminé" >> "${log_folder}/log_STEP1"
    echo "log_1C: processus terminé" >> "${log_folder}/log_STEP1"
    echo "log_1C" > log_1C
    """
} 


//Etape 1D:
process _1D_miniprotsearch {
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
    conda 'bioconda::cd-hit'

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
    conda 'bioconda::mmseqs2'

    input:
    path log_1F

    output:
    path 'log_2'

    script:
    """
    #Etape 2A: Effectuer la recherche MMseqs2
    mkdir -p " ${db_folder_mmseqs_query}"
    python "${project_dir}/scripts/STEP_2/STEP2A_create_db_mmseqs.py" "${orf_pred_file_derep}/ORF_pred_derep.final" ${db_folder_mmseqs_query}

    #log
    echo "log_2A: creation db query mmseqs : terminé" > "${log_folder}/log_STEP2"   
    
    #Etape 2B: Effectuer la recherche MMseqs2
    mkdir -p "${alignment_temp_folder}"
    python "${project_dir}/scripts/STEP_2/STEP2B_blast_mmseqs.py" ${db_folder_mmseqs_query} ${ref_folder} ${alignment_folder} ${alignment_temp_folder}

    #log
    echo "log_2B: alignement general mmseqs : terminé" > "${log_folder}/log_STEP2"    

    #Etape 2C: Extraire les informations
    mkdir -p ${info_folder_IDS}
    mkdir -p ${info_folder_ID}
    mkdir -p ${info_folder_countdup}
    python "${project_dir}/scripts/STEP_2/STEP2C_extract_info.py" ${alignment_folder} ${info_folder}
    python "${project_dir}/scripts/STEP_2/STEP2C_countduplicate.py" ${info_folder_IDS} ${info_folder_countdup}

    #log
    echo "log_2C: extraction des informations : terminé" >> "${log_folder}/log_STEP2"

    #Etape 2D: Extraire les séquences FASTA à partir des identifiants
    mkdir -p ${orf_codants_folder}
    python "${project_dir}/scripts/STEP_2/STEP2D_extractfastafromID.py" ${orf_pred_file_derep} ${info_folder_ID}  ${orf_codants_folder}

    #log
    echo "log_2D: extraction des seqs fasta : terminé" >> "${log_folder}/log_STEP2"
    echo "log_2: processus : terminé" >> "${log_folder}/log_STEP2"
    echo "log_2" > log_2 
    """
}





			// ETAPE 3 : BACT DNA POL remove 

//Etape 3:
process _3_DNApolthrim {
    conda 'bioconda::mmseqs2'

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
    python "${project_dir}/scripts/STEP_3/STEP3A_DNA_pol_init.py" ${bact_DNA_pol} ${ref_folder}/dnapol.prt ${DNA_pol_fasta}
    echo "log_3A: STEP3_DNA_pol_init : terminé" > "${log_folder}/log_STEP3"    
    
    #Etape 3B:Construction DB par mmseqs
    python "${project_dir}/scripts/STEP_3/STEP3B_create_db_mmseqs.py" ${DNA_pol_fasta} ${DNA_pol_fasta_db}
    echo "log_3B: Construction DB par mmseqs : terminé" >> "${log_folder}/log_STEP3"
 
    #Etape 3C: Alignement DB avec les ORF pred  
    python "${project_dir}/scripts/STEP_3/STEP3C_filter_file.py" ${orf_codants_folder} ${DNA_pol_files}
    
    #Etape 3D : blast 
    python "${project_dir}/scripts/STEP_3/STEP3D_blast_mmseqs.py" ${DNA_pol_fasta_db} ${DNA_pol_files} ${DNA_pol_align} ${DNA_pol_align_temp}
    echo "log_3C: Alignement DB avec les ORF pred : terminé" >> "${log_folder}/log_STEP3"

    #Etape 3E: receuil des infos
    python "${project_dir}/scripts/STEP_3/STEP3E_DNA_pol_filt.py" ${DNA_pol_align} ${DNA_pol_align_alnfilt} ${DNA_pol_align_alnfilt_ID}
    echo "log_3E: receuil des infos : terminé" >> "${log_folder}/log_STEP3"

    #Etape 3F: recuperation des fastas
    python "${project_dir}/scripts/STEP_3/STEP3F_extractfastafromID.py" ${orf_pred_file_derep} ${DNA_pol_align_alnfilt_ID} ${DNA_pol_align_seqfilt}
    echo "log_3F: recuperation des fastas : terminé" >> "${log_folder}/log_STEP3"
        
    #Etape 3G: recuperation des fastas
    python "${project_dir}/scripts/STEP_3/STEP3G_dna_pol_filt_aln.py" ${db_folder_mmseqs_ref}/dnapol ${DNA_pol_align_seqfilt} ${DNA_pol_align_final_aln} ${DNA_pol_align_final_aln_temp}
    echo "log_3G: alignement final avec sequences filtree : termine" >> "${log_folder}/log_STEP3"

    #log
    echo "log_3" > log_3 
    """
}


//Etape 4: Polinto-like decision
process _4_polinto_thrim {
    conda 'bioconda::mmseqs2'

    input:
    path log_3

    output:
    path 'log_4'
    
    script:
    """
    # Etape 4A: Init des db 
    python "${project_dir}/scripts/STEP_4/STEP4A_combine_polinto_ref.py" ${ref_folder} ${polinto_fasta}
    echo "log_4: recuperation des fastas : terminé" > "${log_folder}/log_STEP4"
    
    # Etape 4B: Creation des db 
    python "${project_dir}/scripts/STEP_4/STEP4B_create_db.py" ${polinto_fasta} ${polinto_db}
    echo "log_4: creation des db : terminé" >> "${log_folder}/log_STEP4"

    # Etape 4C: rassembler les ORF polinto  
    python "${project_dir}/scripts/STEP_4/STEP4C_collect_orf.py" ${orf_codants_folder} ${polinto_orf_init}
    echo "log_4: rassembler les orf initiaux : terminé" >> "${log_folder}/log_STEP4"

    # Etape 4D: Blast des orf contre les db
    python "${project_dir}/scripts/STEP_4/STEP4D_blast_mmseqs.py" ${polinto_db} ${polinto_orf_init} ${polinto_aln} ${polinto_aln_temp}
    echo "log_4: rassembler les orf initiaux : terminé" >> "${log_folder}/log_STEP4"

    # Etape 4E: Filtration selon les aln 
    python "${project_dir}/scripts/STEP_4/STEP4E_collect_aln_res.py" ${polinto_aln} ${polinto_orf_init} ${polinto_aln_filt} ${polinto_aln_filt_id} ${polinto_orf_final}
    echo "log_4: Filtration selon les aln : terminé" >> "${log_folder}/log_STEP4"

    # Etape 4F: Alignement des orf finaux
    python "${project_dir}/scripts/STEP_4/STEP4F_final_alignment.py" ${polinto_db} ${polinto_orf_final} ${polinto_aln_final} ${polinto_aln_final_temp}                                                                 
    echo "log_4: creation des db : terminé" >> "${log_folder}/log_STEP4"
    
    #log 
    echo "log_4" > log_4 
    """
}

// Extraction des résultats de Miniprot
process _5_final_treatment {
    input:
    path log_4

    output:
    path 'log_5'

    script:
    """
    #Extraction des résultats miniprot
    python "${project_dir}/scripts/STEP_5/STEP5_final_treatment_orf.py" ${orf_codants_final} ${alignment_final} ${alignment_folder} ${orf_codants_folder} ${DNA_pol_align_final_aln} ${DNA_pol_align_seqfilt} ${polinto_orf_final} ${polinto_aln_final}
    echo "log_5: Extraction des résultats généraux terminée" > "${log_folder}/log_STEP5"

    #log
    echo "log_5" > log_5
    """
}