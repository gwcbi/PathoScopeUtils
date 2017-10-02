# PathoScopeUtils
Utility scripts for analyzing PathoScope output


# Prerequisites

+ R
+ Python
+ samtools
+ Picard



# Example workflow


## Step 0: Initial setup

Pathoscope output for each sample should be contained within its own directory.
Files for each sample should include:

+ `pathoid-sam-report.tsv` - Tab-delimited PathoScope report.
+ `updated_outalign.sam` - Updated alignment file in SAM format.




## Step 1: Add taxonomy information to all reports

For each taxonomy report, add columns corresponding to the taxonomic ranks.

#### Output:

+ `<sample>/pathoid-sam-report.withtax.tsv` - PathoScope report with taxonomy for each row

```bash
cat samples.txt | while read s; do
    add_tax_info.R $s/pathoid-sam-report.tsv > $s/pathoid-sam-report.withtax.tsv
done
```

## Step 2: Sort all the updated SAM files

#### Output:

+ `<sample>/updated_outalign.sorted.bam` - Sorted BAM file with all references
+ `<sample>/updated_outalign.sorted.bam.bai` - BAM index

```bash
cat samples.txt | while read s; do
    samtools view -ub $s/updated_outalign.sam | samtools sort > $s/updated_outalign.sorted.bam
    samtools index $s/updated_outalign.sorted.bam
done
```

## Step 3: Get the top OTUs

#### Output:

+ `top_otus.txt` - Top taxonomy IDs for all samples
+ `top_otus.withtax.txt` - Top taxonomy IDs for all samples, including taxonomic name

```bash
python top_otus.py */pathoid-sam-report.tsv > top_otus.txt
python top_otus.py --with_tax */pathoid-sam-report.withtax.tsv > top_otus.withtax.txt
```

## Step 4: Extract top OTUs from BAM files

#### Output:

+ `otu_analysis/<tax_name>/orig/<sample>.bam` - Sorted BAM for each OTU and sample

```bash
python split_otu.py --taxfile top_otus.withtax.txt */updated_outalign.sorted.bam
```

Relocate these files to `otu_analysis/<tax_name>/orig` directory. (The final BAM files will be in the parent directory).

```bash
for d in otu_analysis/*; do
    mkdir -p $d/orig
    mv $d/*.bam $d/orig
    mv $d/*.bam.bai $d/orig
done
```


## Step 5: Fix SAM flags

PathoScope output does not have the correct flags set for secondary alignments. In order to have better read count and coverage estimation, we want to remove secondary alignments. We achieve this here with the following shell pipeline: **1)** Sort by read name (`samtools sort -n`). **2)** Ensure that each read has only one "primary" alignment and the other alignments are secondary (`fix_sam_flags.py`). **3)** Remove the secondary alignments (`samtools view -F 0x100`) and **4)** sort by coordinate order (`samtools sort`).

#### Output:

+ `otu_analysis/<tax_name>/<sample>.bam` - Sorted BAM with primary alignments only for each OTU and sample.


```bash
for f in otu_analysis/*/orig/*.bam; do
    n=$(dirname $(dirname $f))/$(basename $f)
    samtools sort -n $f |\
        samtools view -h |\
        python fix_sam_flags.py |\
        samtools view -ub -F 0x100 |\
        samtools sort \
        > $n
    samtools index $n
    echo $n
done
```

## Step 6: Visualize contig counts


```bash
for d in otu_analysis/*; do
    echo "#-------- $d --------#"
    Rscript analyze_contig_counts.R $d/*.bam
done

for d in otu_analysis/*; do
    echo "#-------- $d --------#"
    Rscript plot_contig_counts.R $d/*.bam
done

for d in otu_analysis/*; do
    echo "#-------- $d --------#"
    Rscript plot_coverage.R $d/*.bam
done
```