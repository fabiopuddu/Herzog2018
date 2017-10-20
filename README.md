# Herzog2018
#### analysis.sh
This shell script is an example of how the commands are sequenced and details the usage of other perl scripts that we provide.
#### gt-filter-lax.pl
This script contains a custom filter for vcftool's vcf-annotate command which will filter on the GQ tag of vcf file lines.
#### vcf_to_gene_list.pl
This script will take a vcf file containing all samples of the experiment and return a list of mutations resulting in defined consequences (as annotated by the Ensembl VEP) and count how many instances of each mutation was observed across samples.  
#### plot_protein.pl 
This script will go through a list of mutations and will create a graphical represenation of their position using a gnuplot script (plot_protein_map.gpl)
