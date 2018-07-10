# A command line tool for running the annotation workflow
# explained below:
# 	We predicted proteins from all genomes using Prodigal, and subsequently identified
#	protein orthologous groups (OGs) using proteinortho v5.16b with default parameters (56). For
#	each OG we chose the longest member as a representative compared these proteins to the
#	EggNOG release 4.5 (57), Pfam release 31 (58), and TigrFam release 15.0 (59) databases for
#	annotation using HMMER3 (60). For EggNOG we downloaded all NOG hmms from the EggNOG
#	website on February 1st, 2018, and ran hmmsearch with an e-value cutoff of 1e-5. For Pfam and
#	TigrFam annotations we used the the noise cutoffs in each HMM as lower bounds for annotationi

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
