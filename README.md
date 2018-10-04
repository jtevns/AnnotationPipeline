# AnnotationPipeline
An annotation pipeline using snakemake that curently supports pfam annotations

To run:
    -Create a directory for annotation
    -Run filtering script on a directory of bins
        -This is required to filter out bins that do not meet the length requirement for running prodigal (20000 bases)
        - Run: filter_and_link_bins.sh  /path/to/bins/to/annotatate
        - the script results in the creation of a directory called Bins in the directory the script is run from with all bins linked in the           directory that meet the length threshold
    -Run Annotation.smk
        -snakemake -s Annotation.smk --cores 
