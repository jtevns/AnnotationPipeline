# A command line tool for running the annotation workflow
# explained below:
#	1. call genes with prodigal
#	2. select orthologous clusters
#	3. select longest member of cluster as
#	       representative
#	4. hmm search eggNOG, pfam,tigrfam
import click
import subprocess
from os.path import basename 
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
# commands callable by user
@click.group()
def main():
    pass

@main.command()
def annotate_wf():
    pass

@main.command()
@click.argument('fasta', nargs=-1)
def prodigal(fasta):
    #prodigal -d out_nuc.fasta -i inFile.fasta -o outputFile.txt -s geneScores.txt
    for inFasta in fasta:
        inBaseName = basename(inFasta).split(".")[0]
        print("calling genes with prodigal on " + inFasta)
        subprocess.call(['prodigal','-a',inBaseName + '_genes.faa','-i',inFasta,'-o', '-q'])

@main.command()
@click.argument('faas',nargs=-1)
@click.argument('project',nargs=1)
def proteinOrtho(faas,project):
        #../software/proteinortho_v5.16b/proteinortho5.pl -project=Sample_99669_DAStool_bin15 -clean Sample_99669_DAStool_bin15_genes.fna combined_lm_DAStool_bin41_genes.fna
    command = ['../software/proteinortho_v5.16b/proteinortho5.pl','-project='+str(project),'-clean']
    for faa in faas:
        command.append(faa)
    #subprocess.call(command)
    print("making orthologous protein cluster files")
    makeOrthClusters(str(project) + ".proteinortho",project)

@main.command()
@click.argument('faas',nargs = -1)
@click.argument('project',nargs = 1)
def selectRepresentatives(faas,project):
    #make a dictoinary of gene sizes
    geneSizes = dict()
    geneSeqs = dict()
    for faa in faas:
        for record in SeqIO.parse(faa,"fasta"):
            geneSizes[faa+":"+record.id] = len(record.seq)
            geneSeqs[faa+":"+record.id] = record.seq
    #select the longest gene from each cluster
    with open(str(project + "_proteinOrthoClusters.tsv")) as clusters:
        with open(project + "_representatives.faa","w+") as outFile:
           for cluster in clusters:
               longestID = "unset"
               longestLen = 0
               #lengths = [geneSizes[gene] for gene in cluster.split()]
               for gene in cluster.split("\t"):
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
                           longestLen = currLen
               SeqIO.write(SeqRecord(geneSeqs[longestID],id = longestID,description=""),outFile,'fasta')
@main.command()
def hmmer():
    pass

# other utility commands 
def makeOrthClusters(proteinClusters,project):
    with open(proteinClusters) as clusters:
        with open(str(project) + "_proteinOrthoClusters.tsv","w+") as outfile:
            for line in clusters:
                if(line.startswith("#")):
                    header = line.split()[4:-1]
                else:
                    clusterLine = line.split("\t")[3:-1]
                    combinedClusterRow = list()
                    for lineLoc in range(0,len(clusterLine)): 
                        clusterRow = header[lineLoc] +":"+ clusterLine[lineLoc]
                        if "*" not in clusterRow:
                            combinedClusterRow.append(clusterRow)
                    outfile.write("\t".join(combinedClusterRow) + "\n")
            
            

if __name__ == '__main__':
    main()
