## Dockerfile
FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty
LABEL Description="tophat-2.1.0.Linux update"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git wget samtools unzip python

# Download the tophat-2.1.0
RUN wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz
RUN tar xvf tophat-2.1.0.Linux_x86_64.tar.gz

# Download the bowtie-2.2.5
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip
RUN unzip bowtie2-2.2.5-linux-x86_64.zip

# Clone the tophat-2.2.1.SE.PE wrapper script
ENV TOPH2GIT https://github.com/upendrak/tophat-2.1.0.git
RUN git clone $TOPH2GIT

# Change the permissions and the path for the wrapper script
RUN chmod +x /tophat-2.1.0/tophat-2.1.0.SE.PE.pl && cp /tophat-2.1.0/tophat-2.1.0.SE.PE.pl /usr/bin 
ENV PATH /tophat-2.1.0.Linux_x86_64/:$PATH
ENV PATH /bowtie2-2.2.5/:$PATH

# Entrypoint
ENTRYPOINT ["/usr/bin/tophat-2.1.0.SE.PE.pl"]
CMD ["-h"]

# Build dockerimage from dockerfile
# sudo docker build -t"=ubuntu/cufflinks-2.2.1:latest" .

# cd to sample_files directory and then do test run
# sudo docker run --rm -v $(pwd):/home/upendra_35/cufflinks-2.2.1/sample_files -w /home/upendra_35/cufflinks-2.2.1/sample_files ubuntu/cufflinks-2.2.1:latest --infile SRR070570_WT.fastq.tophat.bam 

# Push the dockerfile to the gitrepo
