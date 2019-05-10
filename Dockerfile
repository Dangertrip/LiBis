FROM conda/miniconda3-centos7

WORKDIR /data/
ENV URL=https://github.com/Dangertrip/LIBisPipe.git

RUN     cd / && \
        yum install -y wget && \
        yum install -y unzip && \
	wget https://github.com/Dangertrip/LIBisPipe/archive/master.zip && \
        unzip master.zip && \
	conda install -c bioconda samtools=1.1 && \
	conda install -c bioconda bedtools && \
        conda install -c bioconda fastqc && \
        conda install -c bioconda cutadapt && \
        conda install -c bioconda trim-galore && \
	pip install -r LIBisPipe-master/requirements.txt && \
	chmod +x LIBisPipe-master/LIBisPipe && \
    cd /data/ 	
	
ENV PATH /LIBisPipe-master:$PATH


