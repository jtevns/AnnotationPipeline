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
        subprocess.call(['prodigal','-a',inBaseName + '_genes.faa','-i',inFasta,'-o',inBaseName+'_genes.gbk','-s',inBaseName + "_geneInfo.txt", '-q'])

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
def selectRepresentative():
    pass

@main.command()
def hmmer():
    pass

# other utility commands 
def makeOrthClusters(proteinClusters,project):
    with open(proteinClusters) as clusters:
        with open(str(project) + "_proteinOrthoClusters.tsv","w+") as outfile:
            clusterCount = 1
            for line in clusters:
                if(line.startswith("#")):
                    header = line.split()[4:-1]
                else:
                    clusterLine = line.split("\t")[3:-1]
                    for lineLoc in range(0,len(clusterLine)): 
                        clusterRow = str(clusterCount) +"\t"+ header[lineLoc] +"\t"+ clusterLine[lineLoc]
                        if "*" not in clusterRow:
                            outfile.write(clusterRow + "\n")
                clusterCount = clusterCount + 1
            
            

if __name__ == '__main__':
    main()
