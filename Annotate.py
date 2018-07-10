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
        print("calling genes with prodigal")
        subprocess.call(['prodigal','-d',inBaseName + '_genes.fna','-i',inFasta,'-o',inBaseName+'_genes.gbk','-s',inBaseName + "geneInfo.txt", '-q'])

@main.command()
@click.argument('fnas',nargs=-1)
@click.argument('project',nargs=1)
def protienOrtho(fnas,project):
        #../software/proteinortho_v5.16b/proteinortho5.pl -project=Sample_99669_DAStool_bin15 -clean Sample_99669_DAStool_bin15_genes.fna combined_lm_DAStool_bin41_genes.fna
        subprocess.call(['../software/proteinortho_v5.16b/proteinortho5.pl','-project',project,'-clean'])

@main.command()
def selectLongest():
    pass

@main.command()
def hmmer():
    pass

if __name__ == '__main__':
    main()
