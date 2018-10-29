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
        
#function to map to clusters

def MapToClusters(clusters,annotation):
    with open(clusters) as clusters:
        clusterDict = {x.split("\t")[0]:x.split("\t")[1:] for x in clusters }

    outFileName = "MappedToClusters_"+annotation[:-3]+"csv"

    with open(annotation) as annotations, open(outFileName,"w+") as outFile:
        for line in annotations:
            if "#" not in line:
                clusterNum = line.split(" ")[0].split("_cluster_")[1]
                foundClust = clusterDict[clusterNum]
                for geneList in foundClust:
                    if len(geneList.split(",")) > 1:
                        geneList = geneList.split(",")
                        prefix = geneList[0].split(":")[0]
                        genes =  [":".join([prefix,x]) for x in geneList[1:]]
                        for gene in genes:
                            outLine = ",".join([gene.strip()]+line.split()[1:])
                            outFile.write(outLine + "\n")
                    else:
                        outLine = ",".join([geneList.strip()]+line.split()[1:])
                        outFile.write(outLine+"\n")

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

# function for selecting representatives
def selectRepresentatives(faas,clusters):
    #make a dictoinary of gene sizes    
    geneSizes = dict()    
    geneSeqs = dict()    
    for faa in faas:
        for record in SeqIO.parse(faa,"fasta"):
            geneSizes[basename(faa)+":"+record.id] = len(record.seq)
            geneSeqs[basename(faa)+":"+record.id] = record.seq
    #select the longest gene from each cluster
    clusters = open(clusters,"r")
    outFile = open("representatives.faa","w+")
    for cluster in clusters:
        clusterNum = cluster.split("\t")[0]
        longestID = "unset"
        longestLen = 0
        #lengths = [geneSizes[gene] for gene in cluster.split()]
        for gene in cluster.split("\t")[1:]:
            temp = gene.strip().split(":")
            binID = temp[0]
            if "," in temp[1]:
                contigs = temp[1].split(",")
                for contig in contigs:
                    currContig = binID+":"+contig
                    currLen = geneSizes[currContig]
                    if currLen > longestLen:
                        longestID = currContig
                        longestLen = currLen
            else:
                currContig = gene.strip()
                currLen = geneSizes[currContig]
                if currLen > longestLen:
                    longestID = currContig
                    longestIDwNum = currContig +"_cluster_"+str(clusterNum)
                    longestLen = currLen
        SeqIO.write(SeqRecord(geneSeqs[longestID],id = longestIDwNum,description=""),outFile,'fasta')

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

#function to make cluster file
def makeOrthClusters(proteinClusters):
    clusterCount = 0
    clusters = open(proteinClusters, 'r')
    outfile = open("ProteinOrthoClusters.tsv","w+")
    for line in clusters:
        if(line.startswith("#")):
            header = line.split()[4:]
        else:
            clusterLine = line.split("\t")[3:]
            combinedClusterRow = list()
            for lineLoc in range(0,len(clusterLine)):
                clusterRow = header[lineLoc] +":"+ clusterLine[lineLoc]
                if "*" not in clusterRow:
                    combinedClusterRow.append(clusterRow)
            outfile.write(str(clusterCount)+ "\t" + "\t".join(combinedClusterRow) + "\n")
            clusterCount += 1
    clusters.close()
    outfile.close()
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
