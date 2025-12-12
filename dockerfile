FROM continuumio/miniconda3

ENV CACHE_TTL=15
ENV CACHE_MAX_SIZE=1024

EXPOSE 80

COPY ./app /workspace/app
COPY ./environment.yml /workspace/
COPY ./data/eval /workspace/data/eval
COPY ./data/ncbi_lineages_2023-06-15.csv.gz /workspace/data/ncbi_lineages_2023-06-15.csv.gz

WORKDIR /workspace
RUN mkdir /logs

RUN conda env update --file environment.yml
