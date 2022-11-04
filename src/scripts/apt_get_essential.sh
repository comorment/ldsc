#!/bin/bash
apt-get update && apt-get install -y  --no-install-recommends python3=3.8.2-0ubuntu2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

update-alternatives --install /usr/bin/python python /usr/bin/python3 10
