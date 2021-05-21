# AnnotationPipeline
An annotation pipeline using snakemake that curently supports pfam, eggnog, tigrfam annotations. 

To run:  
    -Create a directory for annotation  
    -Create a folder and put your original fasta files in
    - module load python3.7-anaconda
    -Run filtering script on a directory of bins  
 &nbsp;&nbsp;&nbsp;&nbsp;-This is required to filter out bins that do not meet the length requirement for running prodigal (20000 bases and prep the directory for the pipeline)   
 &nbsp;&nbsp;&nbsp;&nbsp;- Run: ```/nfs/turbo/lsa-dudelabs/pipelines/AnnotationPipeline/scripts/prep_annotation.py folder_of_fastas```  
 &nbsp;&nbsp;&nbsp;&nbsp;- the script results in the creation of a directory called Passing_Bins in the directory the script is run from with all 
 &nbsp;&nbsp;&nbsp;&nbsp;bins linked in the directory that meet the length threshold.  
    - copy the config template from /nfs/turbo/lsa-dudelabs/pipelines/AnnotationPipeline/config_template.yaml and name it config.yaml.  
    - use the file extension in the Passing_Bins dir (fa or fna or fasta etc).  
    - module load singularity.  
    -Run Annotation.smk  (I recommend putting this in a slurm script with 36 cores and 180gb mem for fastest run time).  
    &nbsp;&nbsp;&nbsp;&nbsp; ```singularity exec /nfs/turbo/lsa-dudelabs/pipelines/AnnotationPipeline/Singularity/annotation_tools.sif snakemake -s /nfs/turbo/lsa-dudelabs/pipelines/AnnotationPipeline/Annotation.smk --cores```   
Note on headers:
fasta headers cannot contain ":" 
