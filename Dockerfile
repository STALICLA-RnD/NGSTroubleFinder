
# Base image
FROM python:3.10.10-slim-bullseye as base

COPY . /app
WORKDIR /app

## Execute any command that is needed in the Dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends gcc libhts-dev build-essential && rm -rf /var/lib/apt/lists/*

#Install the requirements and the module
RUN python3 -m pip install --upgrade --no-cache-dir pip && \
  python3 -m pip --no-cache-dir install --upgrade build && \
  python3 -m build && \
  python3 -m pip install --no-cache-dir .

## Change following labels as pertinent, this is only informative.
LABEL org.label-schema.maintainer="Samuel Valentini"
LABEL org.label-schema.email="samuel.valentini@stalicla.com"
LABEL org.label-schema.description="ngsTroubleFinder"
# "${MAIN_TAG}"

## When the image is called with "docker run [...]" it will automatically become "docker run [...] ENTRYPOINT", which should be the entry script or binary used
CMD ["ngsTroubleFinder"]

