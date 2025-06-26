# *S. aureus* CC59 genomic analysis #
This repository outlines the analysis pipeline for the paper *S. Jiang. et al. Global emergence and evolution of Staphylococcus aureus clonal complex 59 （in review)*

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly --careful -k 21,33,55,77,99,127
```
Species confirmation by GTDB-Tk (https://github.com/Ecogenomics/GTDBTk)  
```
gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk/classify --cpus 10 --skip_ani_screen
```
QC assessment by checkM (https://github.com/Ecogenomics/CheckM)  
```
checkm lineage_wf -x fasta input_bins output_folder
```
## Genome annotation ##  
Virulence factors (VFs) identification using curated *S. aureus* VF data ```curated_SAVF.faa``` by abriate (https://github.com/tseemann/abricate)  

Antibiotic resistance genes (ARGs) identification using RGI (https://github.com/arpcard/rgi?tab=readme-ov-file)
```
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json /path/to/card.json --local
rgi main --input_sequence nucleotide_input.fasta --output_file output_file --local --clean
```
SCCmec typing and spatyping by bactopia (https://bactopia.github.io/v3.0.0/bactopia-tools/staphopiasccmec/)
```
staphopia-sccmec --assembly test/GCF_001580515.1.fna
bactopia --wf spatyper --bactopia /path/to/your/bactopia/result
```

Pangenome analysis using prokka (https://github.com/tseemann/prokka) and roary (https://sanger-pathogens.github.io/Roary/)  
```
prokka /path/to/"$sample".fasta --quiet --outdir /path/to/prokka_output/"$sample" --force --prefix $sample
roary –f output_dir *.gff
```
## Phylogenetic analysis ##
core genome alignment was generated using snippy (https://github.com/tseemann/snippy), and recombination sites were removed with Gubbins (https://github.com/nickjcroucher/gubbins). A maximum-likelihood phylogenetic tree was then constructed using IQ-TREE (http://www.iqtree.org/) based on clean core genome SNP alignments.
```
snippy --outdir mut1 --ref ref.gbk --ctgs mut1.fasta
run_gubbins.py -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
iqtree -s clean.core.aln --boot-trees --wbtl -m GTR+I+G -B 1000 -nt 18
```
## Genomic population ##
```
library(fastbaps)
library(ape)

fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
baps.hc <- fast_baps(sparse.data)
clusters <- best_baps_partition(sparse.data, as.phylo(baps.hc))
```
## GWAS analysis ##
The GWAS was conducted simultaneously on the pangenome and k-mer association.
Pangenome was retrieved from roary and the gene_presence_absence.csv and traits file were input into Scoary (https://github.com/AdmiralenOla/Scoary).
```
scoary -g gene_presence_absence.csv -t traits.csv -c BH --no_pairwise
```
Pyseer was used for K-mer-based analysis and the tutorial can be accessed in (https://pyseer.readthedocs.io/en/master/tutorial.html#k-mer-association-with-mixed-effects-model)
 
## Prophage analysis ##
The intact prophage was retrieved from PHASTEST using ```PHATEST_API.sh```
```
bash phastest_api.sh --submitjob --inputDir ./00.fnafiles
bash phastest_api.sh --getresults --outDir ./01.phastest_results
```
Pharokka was used to annotate all intact prophages and the pangenome was generated using Roary.
```
pharokka.py -i prophage.fasta -o outdir -d path/to/database_dir -t 16 --fast
roary –i 90 –f output_dir *.gff
```
The network of community pangenome was generated following (https://github.com/JDHarlingLee/GraPPLE). Briefly, the ```gene_presence_absence.csv``` from roary was converted into binary matrix using ```gene_matrix_to_binary.py``` 
```
python gene_matrix_to_binary.py -i gene_presence_absence.csv -o binary_presc_absc.tsv --start_col 15 --delimiter ","
```
Then, the pairwise similarity between genomes from a binary matrix was calculated using ```pw_similarity.py```
```
python pw_similarity.py -i binary_presc_absc.tsv -o example1 -r "isolates" -s "jaccard" -f 0.5
```
## Evolutionary analysis ##
The temporal origins of CC59 genomes were estimated using BactDating (https://github.com/xavierdidelot/BactDating)
```
library(ape)
tree=read.tree('filename.nwk')
year <- scan(input_file, what = numeric(), quiet = TRUE)
res=bactdate(tree,year)
```
The geographical transmission patterns were esimated using TreeTime with mugration model (https://treetime.readthedocs.io/en/latest/tutorials/mugration.html)
```
treetime mugration --tree tree.nwk --states metadata.csv --attribute country
```
The phylodynamic inference of effective population size was conducted using Skygrowth (https://github.com/mrc-ide/skygrowth)
```
require(skygrowth)
require(ape)
require(ggplot2)
tree <- read.tree("your_tree_file.nwk")
mcmcfit <- skygrowth.mcmc( tree)
```
## Other Analysis ###
For network visualisation, we used Cytoscape (https://github.com/cytoscape/cytoscape). The software can also be used for calculating network descriptive statistics. All plots were generated using ggplot2 (https://ggplot2.tidyverse.org/).
