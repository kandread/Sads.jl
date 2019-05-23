FROM julia:latest
MAINTAINER Kostas Andreadis

RUN apt-get update && apt-get -y install bzip2 #build-essential

COPY deps.jl deps.jl

RUN julia deps.jl

COPY scripts/swot.jl swot.jl

RUN julia -e "using Sads; using NCDatasets; using Distributions"

ENTRYPOINT ["julia", "swot.jl"]

CMD []


