# RNAseq_workshop_helpers
A few hole-filling scripts that were useful in preparing our RNAseq workshop

Be warned all of these scripts are currently of a worked-for-me level of testing.

Call each script with `-h` for more usage information.

## agi_finder.py
This script basically just has a regex to get AGIs (Arabidopsis thaliana) gene
identifiers out of prose or tables.

## clean_cogent_output.py
- requires biopython to be installed

This script just checks the output of [Cogent](https://github.com/Magdoll/Cogent) 
and then concatenates all the 'unassigned transcripts' and Cogent output into one file.

## count_blat.py
This script can be used to count up the BLAT or BLAST hits to particular target
sequences. Mostly useful when one has performed cross-species mapping.

- BLAT must be ran with `-out=blast8`
- BLAST must be ran with `-outfmt=6`

## gff3_to_hints_isoseq.py
- requires [dustdas](https://github.com/janinamass/dustdas) to be installed.

Converts the output created by `collapse_isoforms_by_sam.py` from 
[cDNA_cupcake](https://github.com/Magdoll/cDNA_Cupcake) into hints for Augustus.

Note that `collapse_isoforms_by_sam.py` produces gff, and this script parses gff3.
So you will want to convert in between (e.g. with [gffread](https://github.com/gpertea/gffread))

## subset_genome_related.py
- requires samtools (tested on 1.9) to be installed 

This script can ochestrate the filtering of fasta, gff, and bam files to a sub region.
E.g. when you need smaller test-data for a workshop. 


