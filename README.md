# AnnotationPipeline
An annotation pipeline using snakemake that curently supports pfam annotations  

To run:  
    -Create a directory for annotation  
    -Run filtering script on a directory of bins  
 &nbsp;&nbsp;&nbsp;&nbsp;-This is required to filter out bins that do not meet the length requirement for running prodigal (20000 bases)   
 &nbsp;&nbsp;&nbsp;&nbsp;- Run: filter_and_link_bins.sh  /path/to/bins/to/annotatate  
 &nbsp;&nbsp;&nbsp;&nbsp;- the script results in the creation of a directory called Bins in the directory the script is run from with all    &nbsp;&nbsp;&nbsp;&nbsp;bins linked in the directory that meet the length threshold  
    -Run Annotation.smk  
    &nbsp;&nbsp;&nbsp;&nbsp; -snakemake -s Annotation.smk --cores   
