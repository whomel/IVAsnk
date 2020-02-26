# IVAsnk
Consensus of influenza viruses using IVA and snakemake

`src/who_flu.py` : Given multiple paired fastq files, the pipeline runs a quality check (trimmomatic) and a denovo assembly using IVA. Next, blast contigs and extract the closest reference from NCBI DB -> consensus from the mapped reads -> renames/group the results and create the depth plots using circos. 

`src/who_flu_noInd.py` : Same as `who_flu.py` but instead filtering the INDELS from the vcf files before extract the consensus, it deletes all the INDELS from the bam file.  

**Input files:** (input1) Fastq folder with paired reads.

_Usage_:
`snakemake --snakefile who_flu.py -j 48 --config ifq=/Data/seq/raw/illumina/SeqRun10/ out=out_seq/SeqRun10/ -np`

(-j number of cores; ifq fastq folder; -np dry run)

**Output files:** Assemblies and depth plots.

_Dependencies:_ [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [IVA](https://sanger-pathogens.github.io/iva/), [circos](http://circos.ca/), [FluDB database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database)

**[!!]** _Considerations:_ (1) If no flu is detected in the fastq sample, it creates an empty file, to continue with the pipeline.

Pipeline coded using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and make use of [Python3](https://www.python.org/) and [Biopython](https://biopython.org/).


----------------------------

Copyright (C) February/2020. Vijaykrishna Dhanasekaran and and Miguel Grau López. WHO Collaborating Centre for Reference and Research on Influenza. Disease Ecology and Evolutionary Genetics Laboratory Monash University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
