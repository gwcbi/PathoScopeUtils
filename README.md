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

+ `otu_analysis/<tax_name>/<sample>.bam` - Sorted BAM for each OTU and sample

```bash
python split_otu.py --taxfile top_otus.withtax.txt */updated_outalign.sorted.bam
```

## Step 5: Fix SAM flags

Sort by read name (`samtools sort -n`). Ensure that each read has only one "primary" alignment and the other alignments are secondary, (`fix_sam_flags.py`). Remove the secondary alignments (`samtools view -F 0x100`) and sort by coordinate order.

#### Output:

+ `otu_analysis/<tax_name>/fixed/<sample>.bam` - Sorted BAM with primary alignments only for each OTU and sample.


```bash
for f in otu_analysis/*/*.bam; do
    [[ ! -d $(dirname $f)/fixed ]] && mkdir -p $(dirname $f)/fixed
    n=$(dirname $f)/fixed/$(basename $f)
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
