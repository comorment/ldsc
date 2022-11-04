FROM ubuntu:20.04

ENV TZ=Europe
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp
COPY /scripts/apt_get_essential.sh .
RUN bash apt_get_essential.sh
