FROM conda/miniconda3-centos7

WORKDIR /data/
ENV URL=https://github.com/Dangertrip/LiBis.git

RUN     cd / && \
	    conda install -c bioconda libis && \
    cd /data/ 	
	


