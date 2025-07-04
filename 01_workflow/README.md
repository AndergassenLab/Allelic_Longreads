# Allele-specific Long reads
Workflow for allele-specific analysis in long-read data from PacBio

## ISO-seq
Raw reads have to be demultiplexed and then filtered to remove primer and concatemers. To do so, the following commands have to be used:

`$ lima hifi_reads.bam barcoded_primers.fasta fl.bam --isoseq --peek-guess`

`$ isoseq refine fl.bam primers.fasta flnc.bam`

For further details refer to [pbbioconda](https://github.com/PacificBiosciences/pbbioconda) and [isoseq](https://isoseq.how/) documentations.

## Alignment
The flnc reads are aligned using the PacBio wrapper for minimap2:

`$ pbmm2 align ref.fa flnc.bam aligned.bam --sort -j 4 -J 2`

## Allele-specific Analysis
For more details on the WhatsHap tool, see the [documentation](https://whatshap.readthedocs.io/en/latest/index.html).

The VCF file was generated using the known SNPs between the mice strains ([Keane et al., 2011](https://doi.org/10.1038/nature10413)) and formatted according to the WhatsHap [documentation](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-in-vcfs).

`$ whatshap haplotag -o haplotagged.bam --reference reference.fa snps.vcf aligned.bam --ignore-read-groups --ignore-linked-read --skip-missing-contigs --output-haplotag-list haplotypes.tsv`

`$ whatshap split --output-h1 h1.bam --output-h2 h2.bam aligned.bam haplotypes.tsv`

## Quantification
Reads per genes were quantified using htseq for each haplotype file.


