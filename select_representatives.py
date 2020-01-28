from Bio import SeqIO
from Bio import SeqRecord
from os.path import basename
import sys

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
        longestIDwNum = "unset"
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
                        longestIDwNum = currContig +"_cluster_"+str(clusterNum)
            else:
                currContig = gene.strip()
                currLen = geneSizes[currContig]
                if currLen > longestLen:
                    longestID = currContig
                    longestIDwNum = currContig +"_cluster_"+str(clusterNum)
                    longestLen = currLen
        SeqIO.write(SeqRecord(geneSeqs[longestID],id = longestIDwNum,description=""),outFile,'fasta')

selectRepresentatives(sys.argv[1],sys.argv[2])