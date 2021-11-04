#!/bin/bash

version="v1.0.0"
docker run -it --name mentr_it mkoido/mentr:${version}
# docker rm mentr_it
