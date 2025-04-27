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


#GWAS analysis
prokka --outdir output_folder --force --prefix name --addgenes --compliant --cpus 12 input_fasta
roary -e -mafft -p 12 input_gff
scoary -g gene_presence_absence.csv -t traits_file -p 1E-5 -c BH --no_pairwise -o output_folder
macsyfinder --db-type ordered_replicon --sequence-db input_faster --models CONJScan -o output_folder -w 12
isescan.py --seqfile input_fasta --output output_folder --nthread 12
Prophage date retrival from PHASTEST using PHASTEST_API.sh
mash triangle -E -s 5000 -k 13 input_fasta > /path/to/edgelist.tsv

#Phylogenetic analysis
snippy --outdir output_folder --ref ref.db --ctgs input_fasta --force --cpus 12 --ram 12 --prefix name --report
snippy-core --ref ref.gb input_fasta
snippy-clean_full_aln core.full.aln > clean.full.aln
run_gubbins.py -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
iqtree -s snp.aln.phylip --boot-trees --wbtl -m GTR+I+G -B 1000 -nt 64

#Phylogeographic Analysis
BEAST2 v2.7.3 was used for tree generation with setting of HKY model, optimal relaxed clock, coalescent constant population model and a Markov chain Monte Carlo (MCMC).
Tracer was used for quality evaluation with ESS over 200.
select gene alignment format was converted to PAML format in analysis using dna2paml.py
PAML analysis was conducted using setting file codeml.ctl

