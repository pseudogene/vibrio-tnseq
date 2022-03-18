vibrio-tnseq - Tn-seq pipeline for _Vibrio_ sp.

[![Build Status](https://www.travis-ci.com/pseudogene/vibrio-tnseq.svg?branch=master)](https://www.travis-ci.com/pseudogene/vibrio-tnseq) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/e25f356a3e83436bbedc67440dce958c)](https://www.codacy.com/gh/pseudogene/vibrio-tnseq/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=pseudogene/vibrio-tnseq&amp;utm_campaign=Badge_Grade)

# Vibrio Tn-seq

We foster the openness, integrity, and reproducibility of scientific research.

Scripts and tools used to analyse _Vibrio anguillarum_ Tn-seq project. The data were generated using Illumina MiSeq platform.

## Associated publication

> **Essential genes of _Vibrio anguillarum_ and other _Vibrio_ spp. guide the development of new drugs and vaccines**.
> Bekaert M, Goffin N, McMillan S and Desbois A.
>_Front. Microbiol._ 12:755801

[![DOI](https://img.shields.io/badge/DOI-10.3389%2Ffmicb.2021.755801-blue.svg)](https://doi.org/10.3389/fmicb.2021.755801)


## How to use this repository?

This repository hosts both the scripts and tools used by this study and the raw results generated at the time. Feel free to adapt the scripts and tools, but remember to cite their authors!

To look at our scripts and raw results, **browse** through this repository. If you want to reproduce our results, you will need to **clone** this repository, build the docker, and the run all the scripts. If you want to use our data for our own research, **fork** this repository and **cite** the authors.

## Prepare a docker

All required files and tools run in a self-contained [docker](https://www.docker.com/) image.

### Clone the repository

```sh
git clone https://github.com/pseudogene/vibrio-tnseq.git
cd vibrio-tnseq
```

### Create a docker

```sh
docker build --rm=true -t vibrio-tnseq .

# test
docker run -i -t --rm -t vibrio-tnseq /usr/local/bin/run_pipeline.pl -v --gff /databases/vibrio.gff --cgview 1 --png 1 --infolder /data
```

### Start the docker

To import the raw read files and export the results of your analyse you need to link a folder to the docker. It this example your data will be store in `/home/myaccount` (current filesystem) which will be seem as been `/data` from within the docker by using `-v <USERFOLDER>:/data`.

```sh
docker run -i -t --rm -v <absolute_path>:/data -t vibrio-tnseq /bin/bash
```

### Run manually a new analysis

Make sure your raw read file are in `<absolute_path>`. To run manually a new analysis:

```sh
gunzip -c <input.fastq.gz> > <input.fastq>

cutadapt -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTTCAGAGTTCTACAGTCCGACGATCACAC \
  -a TAACAGGTTGGATGATAAGTCCCCGGTCTCTGTCTCTTATACACATCTCCGAGCCCACGAGAC -O 3 \
  -m 10 -M 18 -e 0.15 --times 2 --trimmed-only \
  -o <output.fastq> \
  <input.fastq>

bowtie2 --no-1mm-upfront --end-to-end --very-fast -x /databases/vibrio -U <output.fastq> -S <output.sam>

sam_to_map.pl --sam <output.sam> --cgview 1 --png 1 -v > output.log
```

### Use your own genome and genome annotation

By default de database available is for _Vibrio anguillarum_ NB10. If you used another species or strain you will have to update the database

```
cd /databases/

# get the genome sequence (FASTA format)
wget -cO - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/786/425/GCF_000786425.1_Vibrio_anguillarum_NB10_serovar_O1/GCF_000786425.1_Vibrio_anguillarum_NB10_serovar_O1_genomic.fna.gz > Vibrio_anguillarum_NB10.fa.gz

# get the genome sequence (GFF format)
wget -cO - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/786/425/GCF_000786425.1_Vibrio_anguillarum_NB10_serovar_O1/GCF_000786425.1_Vibrio_anguillarum_NB10_serovar_O1_genomic.gff.gz > Vibrio_anguillarum_NB10.gff.gz

# pre-process the genome sequence for Bowtie2
gunzip /databases/*.gz
bowtie2-build -q -f Vibrio_anguillarum_NB10.fa Vibrio_anguillarum_NB10

#Your "database" : /database/Vibrio_anguillarum_NB10
#Your "gff": /database/Vibrio_anguillarum_NB10.gff
```

```sh
cd /data/
run_pipeline.pl -v --gff /database/Vibrio_anguillarum_NB10.gff --database /database/Vibrio_anguillarum_NB10 --cgview 1 --png 1 --infolder reads
```

## Run a pipeline

Make sure your compressed raw read files are in `<absolute_path>/reads`. To run a new pipeline:

```sh
run_pipeline.pl --infolder reads
```

### Re-run the analysis of this study

You will need to download the dataset from the EBI ENA repository, project [PRJEB39186](https://www.ebi.ac.uk/ena/browser/view/PRJEB39186):

```sh
docker run -i -t --rm -v <absolute_path>:/data -t vibrio-tnseq /bin/bash

mkdir -p reads
cd reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/[...].fastq.gz
[...]
cd ..
run_pipeline.pl -v --cgview 1 --png 1 --infolder reads
exit
```

## Issues

If you have any problems with or questions about the scripts, please contact us through a [GitHub issue](https://github.com/pseudogene/vibrio-tnseq/issues).
Any issue related to the scientific results themselves must be done directly with the authors.

## Contributing

You are invited to contribute new features, fixes, or updates, large or small; we are always thrilled to receive pull requests, and do our best to process them as fast as we can.

## License and distribution

The content of this project itself including the raw data and work are licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/), and the source code presented is licensed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html).
