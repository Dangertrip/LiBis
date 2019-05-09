FROM conda/miniconda3-centos7

WORKDIR /data/
ENV URL=https://github.com/Dangertrip/LIBisPipe.git

RUN yum install -y git && \
	cd / && \
	git clone $URL && \
	conda install --no-update-dependencies -c bioconda samtools=0.1.19 && \
	conda install --no-update-dependencies -c bioconda bedtools && \
	pip install -r LIBisPipe/requirements.txt && \
	chmod +x LIBisPipe/pl && \
    cd /data/ 	
	
ENV PATH /LIBisPipe:$PATH


