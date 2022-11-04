all: ldsc.sif  # additional .sif files may be listed here for corresponding dockerfiles

%.sif: dockerfiles/%/Dockerfile
	docker build -t $* -f dockerfiles/$*/Dockerfile . && \
	scripts/convert_docker_image_to_singularity.sh $* && \
	scripts/move_singularity_file.sh $*

# Example commands executed by this make file
# 
# Run "make ldsc.sif" will trigger this:
# docker build -t ldsc -f dockerfiles/ldsc/Dockerfile . && \  # build docker container
#     scripts/convert_docker_image_to_singularity.sh ldsc && \         # produces ldsc.sif
#     scripts/scripts/move_singularity_file.sh.sh ldsc        # moves ldsc.sif to ../singularity/.
#  
# Run "make all" will build all singularity containers in this project
