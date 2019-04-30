FROM julia:latest
MAINTAINER Kostas Andreadis

RUN apt-get update && apt-get install bzip2

COPY deps.jl deps.jl

RUN julia deps.jl

COPY scripts/swot.jl swot.jl

ENTRYPOINT ["julia", "swot.jl"]

CMD []


