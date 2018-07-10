# A command line tool for running the annotation workflow
# explained below:
#	1. call genes with prodigal
#	2. select orthologous clusters
#	3. select longest member of cluster as
#	       representative
#	4. hmm search eggNOG, pfam,tigrfam
import click

@click.group()
def main():
    pass

@main.command()
def annotate():
    pass
@main.command()
def prodigal():
    pass

@main.command()
def proteinOrtho():
    pass

@main.command()
def selectLongest():
    pass

@main.command()
def hmmer():
    pass

if __name__ == '__main__':
    main()
