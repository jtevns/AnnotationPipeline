#########################################################
# snakemake workflow for the Annotation of a set of bins
# Author: Jacob Evans
# steps: 
#    - determine orthologous clusters
#    - select a representative from each cluster
#    - search database with representative
#    - map back to all contigs in bins
#    - generate gene count tables
##########################################################
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import basename

BINS = [x for x in os.listdir("./Bins")]

# psuedo rule check if the final genetables are present
rule all:
    input:
        #final output files
        "MappedToClusters_pfam_annotation.csv",
        #"MappedToClusters_tigrfam_annotation.csv",
        #"MappedToClusters_eggnog_annotation.csv"
        
# map annotations to memebers of cluster
rule map_to_clusters:
    input:
        # original annotation tsvs
        annotationPF = "pfam_annotation.txt",
        cluster = "ProteinOrthoClusters.tsv"
     #   "tigrfam_annotation.csv",
    #    "eggnog_annotation.csv"
    output:
        # MappedToClusters annotation tsvs
        "MappedToClusters_pfam_annotation.csv",
  #      "MappedToClusters_tigrfam_annotation.csv",
 #       "MappedToClusters_eggnog_annotation.csv"
    run:
        MapToClusters(input.cluster,input.annotationPF)

# search databases
rule search_databases:
    input:
        #representative fasta
        "representatives.faa"
    output:
        # original annotation tsvs
        "pfam_annotation.txt",
  #      "tigrfam_annotation.csv",
 #       "eggnog_annotation.csv"
    threads:40
    shell:
        "hmmsearch -o pfamlog --cpu {threads} --cut_nc --tblout pfam_annotation.txt /nfs/vdenef-lab/Shared/Jacob/Databases/PFAM-A_hmm/Pfam-A.hmm {input}"


#select representatives
rule select_representatives:
    input:
         #cluster tsv
        clus = "ProteinOrthoClusters.tsv",
        faas = expand("Gene_Calls/{binName}.faa",binName=BINS)
    output:
        #representative fasta
        "representatives.faa"
    run:
        selectRepresentatives(input.faas,input.clus)

# make cluster file from ortho
rule gen_cluster_file:
     input:
         "ProteinOrthoOut.proteinortho"
     output: 
         "ProteinOrthoClusters.tsv"
     run:
         makeOrthClusters(input[0])

# cluster genes from all gene calls
rule cluster_genes:
    input:
        expand("Gene_Calls/{binName}.faa",binName=BINS)
    output:
        "ProteinOrthoOut.proteinortho"
    shell:
        """
        perl /nfs/vdenef-lab/Shared/Jacob/software/proteinortho_v5.16b/proteinortho5.pl -project=ProteinOrthoOut -clean {input}
        """
#call genes from all bins
rule call_genes:
    input:
        "Bins/{binName}"
    output:
        faa="Gene_Calls/{binName}.faa",
        gbk="Gene_Calls/{binName}.gbk"
    shell:
        " prodigal -a {output.faa} -i {input} -o {output.gbk} -q "

#filtering step 
