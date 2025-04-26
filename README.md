# S. aureus CC59 genomic analysis #
This repository outlines the analysis pipeline for the paper S. Jiang. et al. Global emergence and evolution of Staphylococcus aureus clonal complex 59 （in review)

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly --careful -k 21,33,55,77,99,127
```
Species confirmation by GTDB-Tk (https://github.com/Ecogenomics/GTDBTk)  
```
!\gtdbtk classify --genome_dir genomes --align_dir /tmp/gtdbtk/align --out_dir /tmp/gtdbtk/classify -x gz --cpus 2 (--skip_ani_screen
```



#ARGs and VF identification
abricate --db virulence_factor_db input_fasta > vf.tab
resfinder.py -o output_folder -s "Staphylococcus aureus" --acuired --point -ifa input_fasta
The output results were combined and converted using arg_trans.py

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

