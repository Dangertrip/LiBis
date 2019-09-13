FROM conda/miniconda3-centos7

WORKDIR /data/
ENV URL=https://github.com/Dangertrip/LiBis.git

RUN     cd / && \
        yum install -y wget && \
        yum install -y unzip && \
	wget https://github.com/Dangertrip/LiBis/archive/master.zip && \
        unzip master.zip && \
	conda install -c bioconda samtools=1.1 && \
	conda install -c bioconda bedtools && \
        conda install -c bioconda fastqc && \
        conda install -c bioconda cutadapt && \
        conda install -c bioconda trim-galore && \
	pip install -r LiBis-master/requirements.txt && \
	chmod +x LiBis-master/LiBis && \
    cd /data/ 	
	
ENV PATH /LiBis-master:$PATH


