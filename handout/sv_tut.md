Key Learning Outcomes
---------------------

By the end of the structural variant (SV) detection practical course
participants will:

-   Have been provided with key fundamentals on how paired-end mappings
    and split-read/soft-clipped read patterns are used in detecting
    deletions, tandem duplicates, inversions and translocations.

-   Know what important quality control checks need to be evaluated
    prior to structural variant calling.

-   Have run DELLY on a subset of whole genome next generation
    sequencing data pertaining to a single human tumour with a matched
    normal control.

-   Be able to filter high confidence SV predictions.

-   Have gained basic knowledge to interpret the VCF output provided by
    DELLY.

-   Have used their understanding of distinct SV paired-end mapping and
    soft-clipped read patterns to visually verify DELLY predicted SVs
    using IGV.

Resources You’ll be Using
-------------------------

### Tools Used

DELLY
:   \
    <https://github.com/tobiasrausch/delly>

Samtools
:   \
    <http://sourceforge.net/projects/samtools/files/samtools/1.2>

Tabix
:   \
    <http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2>

Vcftools
:   \
    <https://vcftools.github.io/index.html>

Picard
:   \
    <https://github.com/broadinstitute/picard>

Python2.7.10
:   \
    <https://www.python.org/downloads/release/python-2710/>

PyVCF Banyan numpy
:   \
    <https://pypi.python.org/pypi>

Useful Links
------------

SAM Specification
:   \
    <http://samtools.sourceforge.net/SAM1.pdf>

Explain SAM Flags
:   \
    <http://picard.sourceforge.net/explain-flags.html>

Alignment Quality Control
-------------------------

For structural variant calling several alignment quality control metrics
should be evaluated. All paired-end mapping methods heavily rely on the
insert size distribution. GC-content biases are important as it can
impact read-depths. DELLY generates read-depth ratios between tumour and
control samples. The percentage of mapped reads, singletons, duplicates
and properly paired reads are additional metrics you should evaluate
prior to any structural variant calling. These statistics vary largely
by protocol and hence, it is usually best to compare multiple different
sequencing runs using the same against each other to highlight outliers.

It is recommended that Picard module commands `CollectInsertSizeMetrics`
and `CollectGcBiasMetrics`, and `samtools flagstat` command be used.

Prepare the Environment
-----------------------

As a quick introduction we will do a structural variant analysis using a
single immortal cancer cell line and its control genome (mortal parental
cells). Total execution time to run the DELLY structural discovery
calling program for a matched normal tumour pair will vary depending on
the read coverage and the size of the genome. As a guide, it can take
approximately 10 to 50 hours (translocation predictions taking the
longest), for a matched normal tumour pair (each 40-50x coverage)
running on 2 cpus on a server with sufficient RAM.

The bam files we will be working on are a subset of the original WGS bam
files, limited to specific chromosomal regions to speed up the analysis
and to meet the time constraints for this practical.

Firstly, we will use shell variables to help improve the readability of
commands and streamline scripting. Each distinct variable will store a
directory path to either, the input WGS bam files, hg19 reference,
programs or output.

Open the Terminal.

First, go to the right folder, where the data are stored.

    cd /home/trainee/sv
    ls
    mkdir <YourFirstName>
    cd <YourFirstName>

    export DS=/home/trainee/sv/data
    export RF=/home/trainee/sv/reference_data
    export SF=/home/trainee/sv/variantFiltering/somaticVariants
    export BR=/home/trainee/snv/Applications/igv
    export CZ=/home/trainee/sv/converter

Somatic Structural Variant Discovery using DELLY
------------------------------------------------

In order to generate *putative* somatic SVs it is crucial to account for
germline SVs. To facilitate, DELLY requires the joint input of a match
normal control and the cancer aligned sequencing data (bam files).

    delly -t DEL -x $RF/hg19.excl -o del.vcf -g $RF/hg19.fa $DS/cancer_cell_line.bam $DS/control.bam
    delly -t DUP -x $RF/hg19.excl -o dup.vcf -g $RF/hg19.fa $DS/cancer_cell_line.bam $DS/control.bam
    delly -t INV -x $RF/hg19.excl -o inv.vcf -g $RF/hg19.fa $DS/cancer_cell_line.bam $DS/control.bam
    delly -t TRA -x $RF/hg19.excl -o tra.vcf -g $RF/hg19.fa $DS/cancer_cell_line.bam $DS/control.bam

Description of the arguments used in the command:

DEL
:   = conduct deletion discovery

DUP
:   = conduct tandem duplication discovery

INV
:   = conduct inversion discovery

TRA
:   = conduct translocation discovery

-o
:   = vcf output

-g
:   = reference genome in FASTA format

-x
:   = genomic regions to exclude (e.g. centro- and telomeric regions)

DELLY VCF output
----------------

A VCF file has multiple header lines starting with the hash \# sign.
Below the header lines is one record for each structural variant. The
record format is described in the below table:

[H]

lll

Column & Field & Description\
1 & CHROM & Chromosome name\
2 & POS & 1-based position. For an indel, this is the position preceding
the indel\
3 & ID & Variant identifier\
4 & REF & Reference sequence at POS involved in the variant\
5 & ALT & Comma delimited list of alternative sequence(s)\
6 & QUAL & Phred-scaled probability of all samples being homozygous
reference\
7 & FILTER & Semicolon delimited list of filters that the variant fails
to pass\
8 & INFO & Semicolon delimited list of variant information\
9 & FORMAT & Colon delimited list of the format of sample genotypes in
subsequent fields\
10+ & & Individual genotype information defined by FORMAT\

[tab:vcfoutput]

You can look at the header of the vcf file using grep, `-A 1` includes
the first structural variant record in the file:

    grep "^#" -A 1 del.vcf

The INFO field holds structural variant site information whereas all
genotype information (annotated as per the FORMAT fields) is provided in
the sample column. Reference supporting reads are compared to
alternative supporting reads and mapping qualities are used to compute
genotype likelihoods (GL) for homozygous reference (0/0), heterozygous
reference (0/1) and homozygous alternate (1/1) (GT). The final genotype
(GT) is simply derived from the best GL and GQ is a phred-scaled
genotype quality reflecting the confidence in this genotype. If GQ\<15
the genotype is flagged as LowQual. The genotyping takes into account
all paired-ends with a mapping quality greater than 20 by default.

The INFO field provides information on the quality of the SV prediction
and breakpoints. If you browse through the vcf file you will notice that
a subset of the DELLY structural variant predictions have been refined
using split-reads. These precise variants are flagged in the vcf info
field with the tag *PRECISE*. To count the number of precise and
imprecise variants you can simply use grep.

    grep -c -w "PRECISE" *.vcf
    grep -c -w "IMPRECISE" *.vcf

DELLY clusters abnormal paired-ends and every single cluster gives rise
to an IMPRECISE SV call. For every IMPRECISE SV call an attempt is made
to identify supporting split-reads/soft-clipped reads DELLY then
computes a consensus sequence (INFO:CONSENSUS) out of all split-read
candidates and then aligns this consensus sequence to the reference
requiring at least -m many aligned bases to the left and right (default
is 13). INFO:PE is the number of supporting paired-ends. INFO:CT refers
connection types (CT), which indicates the order and orientation of
paired-end cluster mappings (e.g. 3to3 for 3’ to 3’ and 5to5 for 5’ to
5’). Values can be 3to5, 5to3, 3to3 or 5to5. Different names exist for
these connection types in the literature, head-to-head inversions,
tail-to-tail inversions, and so on. The consensus alignment quality
(SRQ) is a score between 0 and 1, where 1 indicates 100% identity to the
reference. Nearby SNPs, InDels and micro-insertions at the breakpoint
can lower this score but only for mis-assemblies it should be very poor.
DELLY currently drops consensus alignments with a score\<0.8 and then
falls back to an IMPRECISE prediction.

SVs are flagged as FILTER:LowQual if PE\<3 OR MAPQ\<20 (for
translocations: PE\<5 OR MAPQ\<20), otherwise, the SV results in a
FILTER:PASS. PRECISE variants will have split-read support (SR\>0).

Somatic Structural Variant Filtering
------------------------------------

Please note that this vcf file contains germline and somatic structural
variants but also false positives caused by repeat induced mis-mappings
or incomplete reference sequences. As a final step we have to use the
structural variant site information and the cancer and normal genotype
information to filter a set of confident somatic structural variants.
DELLY ships with a somatic filtering python script. For a set of
confident somatic calls one could exclude all structural variants
\<400bp, require a minimum variant allele frequency of 10%, no support
in the matched normal and an overall confident structural variant site
prediction with the VCF filter field being equal to PASS.

    python $SF/somaticFilter.py -t DEL  -T cancer_cell_line -N control -v del.vcf -o del.filt.vcf -a 0.1 -m 400 -f
    python $SF/somaticFilter.py -t DUP -T cancer_cell_line -N control -v dup.vcf -o dup.filt.vcf -a 0.1 -m 400 -f
    python $SF/somaticFilter.py -t INV  -T cancer_cell_line -N control -v inv.vcf -o inv.filt.vcf -a 0.1 -m 400 -f
    python $SF/somaticFilter.py -t TRA -T cancer_cell_line -N control -v tra.vcf -o tra.filt.vcf -a 0.1 -m 400 -f

Using VCFtools we can merge all somatic structural variants together in
a single vcf file.

    vcf-concat del.filt.vcf dup.filt.vcf inv.filt.vcf tra.filt.vcf | vcf-sort > somatic.sv.vcf

For large VCF files you should also zip and index them using bgzip and
tabix. Please run the below commands to meet the requirements for
visualising somatic structural variants using IGV.

    bgzip somatic.sv.vcf
    tabix somatic.sv.vcf.gz

Visualisation of Somatic Structural Variants
--------------------------------------------

The final step will be to browse some of these somatic structural
variants in IGV and to visually verify the reliability of the calls. VCF
file structure was designed by the 1000Genomes Project consortium with
considerable focus being placed on single-nucleotide variants and
InDels. It is arguable whether the file format is easy to interpret for
basic SVs and much less for complex SV formatting. To make it easier to
see the breakpoints we will create a bed file.

    $CZ/sv.vcf2bed.sh somatic.sv.vcf.gz > somatic.sv.bed
    head somatic.sv.bed

Load the IGV browser

    $BR/igv.sh

Once IGV has started use *File* and *Load from File* to load the
cancer\_cell\_line.bam. To dedicate more of the screen to read pile ups,
right click on the left hand box that contains the
*cancer\_cell\_line.bam Junctions* then select *Remove Track*. Go to
*View* menu and select *Preferences*, then click on the Alignments tab
and in the *Visibility range threshold (kb)* text box, enter 600. This
will allow you to increase your visibility of pile ups as you zoom out.
Now look for check box for *Filter secondary alignments* click on to
ensure we do not see secondary alignments (alternate mapped position of
a read) . Also ensure that *Show soft-clipped bases* has been checked
then click *OK*.

Then import `somatic.sv.bed` from your working directory using *Regions*
and *Import Regions*.

### Verify Deletion

The somatic structural variants can then be browsed easily using the
*Region Navigator*. Select the deletion (chrX:76853017-77014863) from
the *Region Navigator* and click *View*. This will centre the IGV
alignment view on the selected structural variant. Close the regions of
interest pop up window. The red bar below the ruler marks the region of
the deletion.

It’s usually best to zoom out once by clicking on the ’-’ sign in the
toolbar at the top, to give a wide view of the supporting abnormal
paired-end read mappings.

To highlight the abnormal paired-ends right click on the main track
display panel and select *Color alignments by* and then switch to
*insert size and pair orientation*.

Read pairs that have a larger than expected insert size will be
highlighted in red. Click *View as pairs*. Right click and *Sort
alignments by* then select *start location*.

How many abnormal *paired-end read pairs* (red coloured F/R oriented
read pairs) can you see that spans the deletion region? Does this number
coincide with the INFO:PE?

19, YES

Zoom into left breakpoint and tally the number of soft-clipped reads at
both junctions. How many abnormal split-reads (soft-clipped reads) did
you observe? Clue INFO:SR.

11

Go to the RefSeq genes track at the bottom of IGV and right click to
*Expanded*. Is the predicted deletion likely to have a deleterious
impact on a gene? If so, what gene and exons are deleted?

ATRX, exons 2 to 25.

Does this region appear to be completely removed from this cancer
genome? How can you tell?

Yes, there is no read coverage within this deletion region relative to
the control genome.

### Verify translocation

Select the duplication (chrX:45649874-45689322) from the *Region
Navigator*.

Highlight the abnormal paired-ends by clicking and selecting *Color
alignments by* and then switch to *insert size and pair orientation*.
Also, invoke *Sort alignments by* then select *start location*.

Zoom out until you can see all the purple reads at the junction.

What is the direction of the purple cluster of reads (indicates that
mate reads are mapped to Chr15)? Is it pointing to the tail or head of
Chr18?

Forward, towards the tail, or 3’ (+ive)

Right click on to one of the purple coloured reads and select *View mate
region in split-screen*. This will split the screen and display Chr15 on
the left and place a red highlighted outline on both reads, to indicate
the pairs. Right click and *Sort alignments by* then select *start
location*. What is the direction of the yellow cluster of reads
(indicates that mate reads are mapped to Chr18)? Is it pointing to the
tail or head of Chr15? If you wish to zoom in or out, first click inside
of the chromosome ideogram panel, then ctrl- to zoom out and shift+ to
zoom in.

? Reverse, towards the head, or 5’ (+ive)

How is the Chr15 and Chr18 fused (which one of the four translocation
connection types)? If you are uncertain then run a BLAT search using the
INFO:CONSENSUS sequence (gedit somatic.sv.vcf).
https://genome.ucsc.edu/cgi-bin/hgBlat?command=start

RF, head to tail, or 5 to 3. Therefore, Chr18 left side is fused to
Chr15 right side.

Did DELLY predict a reciprocal translocation? How can you tell?

No, as we would expect to observe a Chr18 right side fused to Chr15 left
side, near the same breakpoints.

What gene structures is this translocation predicted to impact?

PARD6G on Chr15 and ADAMTSL3 on Chr18.

What is one possible reason why there is no observable read coverage
after Chr18 breakpoint?

Chromosome loss.

### Verify tandem duplication

Select the tandem duplication (chrX:45649874-45689322) from the *Region
Navigator*.

Highlight the abnormal paired-ends by clicking and selecting *Color
alignments by* and then switch to *insert size and pair orientation*.
Also, invoke *Sort alignments by* then select *start location*.

Zoom out until you can see all the red paired-end reads spanning the two
junctions. After that zoom in on the cluster of abnormal reads on the
left junction and then right junction.

Which is the order and orientation of these paired-end reads (FR, RF, FF
or RR)?

RF

What is the estimated read-depth ratio of the cancer\_cell\_line versus
normal control (INFO:RDRATIO) over the duplicate region?

3.3 ( 3 x increased read depth)

### Verify Inversion

Type into the search box near at tool bar, Chr20:54834492

Right click on main display and select *Group alignments by* then switch
on *paired-orientation*. Also, right click and *Sort alignments by* then
select *start location*.

Zoom out until you can see all the red coloured cluster of reads near
the breakpoint.

Right click on to one of the red coloured reads and select *View mate
region in split-screen*. This will split the screen and display the read
mate on the right side. This will take you to the mate-reads near the
second breakpoint.

Which direction are the paired-end reads spanning (left or right
spanning)?

Right

What is the estimated size of the inverted interval?

55,408,660 – 54,834,492 = 574,168 bp

Acknowledgements
----------------

We would like to thank and acknowledge Tobias Rausch (EMBL Heidelberg)
for his help and for allowing us to borrow and adapt his replies to
questions and original course material.
