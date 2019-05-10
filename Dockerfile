FROM conda/miniconda3-centos7

WORKDIR /data/
ENV URL=https://github.com/Dangertrip/LIBisPipe.git

RUN yum install -y git && \
	cd / && \
	git clone $URL && \
	conda install -c bioconda samtools=0.1.19 && \
	conda install -c bioconda bedtools && \
        conda install -c bioconda fastqc && \
        conda install -c bioconda cutadapt && \
        conda install -c bioconda trim-galore && \
	pip install -r LIBisPipe/requirements.txt && \
	chmod +x LIBisPipe/LIBisPipe && \
    cd /data/ 	
	
ENV PATH /LIBisPipe:$PATH


