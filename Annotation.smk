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
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import basename 


WANTED_ANNOTATIONS = config["ANNOTATIONS"]
BINS = [x for x in glob.glob("Passing_bins/*")]
BIN_NAMES = [x.split("/")[-1].split(".")[0] for x in BINS]

#rule All to indicate ending point
rule all:
    input:
       ["MappedToClusters_{annotation}_annotation.csv".format(annotation=annotation) for annotation in WANTED_ANNOTATIONS]

#call genes from all bins with prodigal
rule call_genes:
    input:
        "Passing_bins/{binName}." + config["INPUT_EXT"]
    output:
        faa="Gene_Calls/{binName}.faa",
        gbk="Gene_Calls/{binName}.gbk"
    shell:
        "prodigal -a {output.faa} -i {input} -o {output.gbk} -q "

# cluster genes from all gene calls with proteinortho
rule cluster_genes:
    input:
        expand("Gene_Calls/{binName}.faa",binName=BIN_NAMES)
    output:
        "annotation.proteinortho.tsv"
    shell:
        """
        proteinortho -project="annotation" -clean {input}
        #mv Protein_Ortho_Out.* Protein_Ortho_Out
        """

# make cluster file from proteinortho output with make_cluster_file.py
rule gen_cluster_file:
     input:
         "annotation.proteinortho.tsv"
     output: 
         "ProteinOrthoClusters_reformatted.tsv"
     script:
         "scripts/make_cluster_file.py"

#select representatives
rule select_representatives:
    input:
        clus = "ProteinOrthoClusters_reformatted.tsv",
        faas = expand("Gene_Calls/{binName}.faa",binName=BIN_NAMES)
    output:
        "representatives.faa"
    script:
        "scripts/select_representatives.py"

rule search_pfam:
    input:
        reps = "representatives.faa",
        db = config["DATABASE_DIRS"]["pfam"]
    output:
        "pfam_annotation.txt"
    threads:36
    shell:
        "hmmsearch -o pfam.log --cpu {threads} --cut_nc --tblout pfam_annotation.txt {input.db} {input.reps}"

rule search_tigrfam:
    input:
        reps = "representatives.faa",
        db = config["DATABASE_DIRS"]["pfam"]
    output:
        "tigrfam_annotation.txt"
    threads:36
    shell:
        "hmmsearch -o tigrfam.log --cpu {threads} --cut_nc --tblout tigrfam_annotation.txt {input.db} {input.reps}"

rule search_eggnog:
    input:
        reps = "representatives.faa",
        db = config["DATABASE_DIRS"]["pfam"]
    output:
        "eggnog_annotation.txt"
    threads:36
    shell:
        "hmmsearch -o eggnog.log --cpu {threads} --cut_nc --tblout eggnog_annotation.txt {input.db} {input.reps}"

# map annotations to members of cluster
rule map_to_clusters:
    input:
        cluster = "ProteinOrthoClusters_reformatted.tsv",
        file = "{annotation}_annotation.txt"
    output:
        "MappedToClusters_{annotation}_annotation.csv"
    script:
        "scripts/map_to_clusters.py"
