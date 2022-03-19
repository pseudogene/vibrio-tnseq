# Copyright 2014-2021, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# This file is part of Vibrio-TNseq.
#
# Vibrio-TNseq is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Vibrio-TNseq is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v3
# along with Vibrio-TNseq. If not, see <http://www.gnu.org/licenses/>.
#
FROM ubuntu:20.04
LABEL description="Vibrio-TNseq Docker" version="1.2" Vendor="Institute of Aquaculture, University of Stirling"

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y python3 python3-pip bowtie2 wget unzip ca-certificates-java default-jre default-jre-headless --no-install-recommends && \
    DEBIAN_FRONTEND=noninteractive apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install cutadapt==3.7

RUN wget -q --no-check-certificate https://github.com/phac-nml/gview-wiki/raw/master/resources/downloads/gview-1.7.zip -O gview-1.7.zip && \
    unzip -q gview-1.7.zip && \
    mv gview/gview.jar /usr/local/bin && \
    rm -rf gview gview-1.7.zip

RUN mkdir /databases
COPY docker/vibrio.fa.gz /databases/vibrio.fa.gz
COPY docker/vibrio.gff.gz /databases/vibrio.gff.gz
COPY docker/sam_to_map.pl /usr/local/bin/sam_to_map.pl
COPY docker/run_pipeline.pl /usr/local/bin/run_pipeline.pl
RUN gunzip /databases/*.gz
RUN bowtie2-build -q -f /databases/vibrio.fa /databases/vibrio

RUN mkdir /data
COPY docker/test.fastq.gz /data/test.fastq.gz

WORKDIR /data
