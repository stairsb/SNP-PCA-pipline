# PCA-pipline
PCA pipeline for whole genome data
read me, how I did it, scripts, figures, metadata


## Table of contents
Obtaining data
Aligning to reference genome


## Obtaining the data
Most of the raw reads were obtained using illumina sequencing data downloaded from the NCBI SRA database. Reads in this analysis that were not downloaded from a database were prepped by the Uehling lab at OSU and sequenced by the UO sequencing core using illumina. The reference genome used was obtained from the NCBI SRA database https://www.ncbi.nlm.nih.gov/Traces/wgs/SMRR01?display=contigs&page=1. 


## Aligning short reads to reference genome
When assembling raw reads there are two options, a denovo assembly or the better method, a reference guided assembly which requires a previously assembled reference genome. Rhizopus microsporus has a well documented reference genome (cited above) so the later of the two were used in this analysis. The refernce genome first has to be indexed using `bwa index` for efficient acess to the short reads within the file. To align reads to a reference genome `bwa mem` was used (BWA documentation http://bio-bwa.sourceforge.net/bwa.shtml).

```
bwa index -a bwstw $REFERENCEGENOME.FA
bwa mem -t 10 $INDEXEDREFERENCEGENOME.FA ${forwardreadsfile1.fq} ${reversereadsfile2.fq} > alignedreads.sam
```

## Organizing the data
We now have SAM files containing the short reads that were aligned. In the current form the downstream analysis will be slow since we using whole genome data for R. microsporus which results in large SAM files. As more genomes are added to the pipline, the time it will take for these programs in our downstrean analysis to run will increase. To combat this we will use `samtools view` to convert our human readable SAM file to binary files which can be process faster and `samtools sort` which will sort alignments by the leftmost coordinates (SAMTOOLS documentation http://samtools.sourceforge.net/samtools-c.shtml). 

```
samtools view -bS alignedreads.sam > alignedreads.bam
samtools sort -@ alignedreads.bam -o sorted.bam
```

## SNP calling
Before moving forward we must create a population map file that contains a list of samples and the sub-populations that each of these samples correspond to. `Gstacks` will find SNPs at each loci for each of the sub-populations and then genotype each sample at the identified SNP. Finally `gstacks` will determine haplotypes using SNPs that were found. (gstacks documentation http://catchenlab.life.illinois.edu/stacks/comp/gstacks.php)

`gstacks -t 8 -I /sorted/bams/directory/ -M pop_map.txt -O /sorted/bams/directory/`

Output:
```
catalog.calls       
catalog.fa.gz 
gstacks.log               
gstacks.log.distribs
```

View output:
```
zcat catalog.calls | less -S
zcat catalog.fa.gz  | less -S
less gstacks.log.distribs
less gstacks.log
```

## Filtering
At this point we have the SNP calls and haplotypes for each sample which are correspended to sub-populations, this would be anough information to input into a PCA analysis but the result may not reflect what actual evolutionary history of our samples are. Depending on the size of our genomes and the amount of samples included here there could future memory and speeds issues when calculating our PCA plot. Adding in data filtering steps will help us keep our analysis accurate and help us eleminate data that we shouldn't actually include. The first program we will use to do this is `populations`. 

Populations filtering:
```
-p — minimum number of populations a locus must be present in to process a locus.
-r — minimum percentage of individuals in a population required to process a locus for that population.
--write-single-snp — restrict data analysis to only the first SNP per locus.
--plink — output genotypes in PLINK format.
```

`populations -t 8 -P /sorted/bams/directory/ -M pop_map.txt -r 0.70 -p 2 --plink --write-single-snp` 

Output:
```
populations.haplotypes.tsv
populations.hapstats.tsv
populations.log
populations.log.distribs
populations.markers.tsv
populations.plink.map
populations.plink.ped
populations.sumstats_summary.tsv
populations.sumstats.tsv
```

The second program we will be using for filtering is `plink` http://zzz.bwh.harvard.edu/plink/thresh.shtml. 'Plink' is a great program for analysis or filtering genotype by sequencing data. 

Plink filtering:
```
--maf - Exclude variants with minor allele frequency lower than a threshold (default 0.01).
--allow-extra-chr - Permit unrecognized chromosome codes.  The '0'
--file - Specify .ped + .map filename prefix (default 'plink').
--geno - Exclude variants with missing call frequencies greater than a threshold (default 0.1).  (Note that the default threshold is only applied if --geno is invoked without a parameter; when --geno is not invoked, no per-variant missing call frequency ceiling is enforced at all.  Other inclusion/exclusion default thresholds work the same way.)
--recode A, count alleles named in the file.
--out - Specify prefix for output files.
```

plink will read both the `.map` and `.ped` files if you input `.plink`. Make sure to give it the directory that both of these files are located in and ask the program to open a `.plink` file

`plink --allow-extra-chr --file populations.plink --maf 0.03 --geno 0.15 --recodeA --out filtered`

Output:
```
filtered.raw
filtered.nosex
filtered.log
```

## PCA analysis in R
The following code will save plots in your current R directory. Some versions of R-Studio will not be able to save the plots this way, you may have to run the following code directly in R to get the output PDFs through this method. You can create a png or jpeg insteads by replacing the code that says "pdf" by your prefered file type.

`library("adegenet")`

Store the data as a genlight object:

`filtered_data <-read.PLINK(file='filtered.raw')`

PCA for large SNPs data sets that are stored as a genlight object:

`SNP_PCA <- glPca(filtered_data, nf = 2)`

Plot will be stored as pdf in current directory:

`pdf(file = "Eigenvalues.pdf")`

Barplot showing eigenvalues:

`barplot(SNP_PCA$eig, main="eigenvalues", col=heat.colors(length(SNP_PCA$eig)))`

Closes the plot:

`dev.off()`

![Eigenvalues](https://user-images.githubusercontent.com/111078377/206596020-8bdea38a-10d6-4588-9962-0d50a607b693.png)


`library("ggplot2")`

Formating the final data before plotting:

```
test <- SNP_PCA$scores |> data.frame()
test2 <- pop(filtered_data)
test$ecology <- test2
```

Store plot as pdf in current directory:

`pdf(file = "PCA_plot.pdf")`

Plotting the data:

`ggplot(test, aes(x=PC1, y=PC2, fill = ecology)) +
  geom_point(size = 3, shape = 21)`

Closes the plot:

`dev.off()`

![PCA_plot](https://user-images.githubusercontent.com/111078377/206596333-1be82136-baa8-4db8-a929-153d243559b0.png)

