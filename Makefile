# Settings
CONDA_ENV=epianeufinder
SHELL=bash
MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh


UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Linux)
    MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    BASE := ${HOME}/miniconda3/bin/conda
endif
ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),x86_64)
        MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    endif
    ifeq ($(UNAME_M),arm64)
        MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-MacOSX-arm64.sh
    endif
endif

ifeq ($(OS),Windows_NT)
    CONDA := $(strip $(shell where.exe conda))
else
    CONDA := $(strip $(shell which conda))
endif


BASE := $(shell dirname $(shell dirname ${CONDA}))

ACTIVATE=${BASE}/bin/activate
DEACTIVATE=${BASE}/bin/deactivate


default: help

install-conda: ## install Miniconda
	echo "installing conda"
	echo base $(BASE)
	echo conda $(CONDA)
	echo activate $(ACTIVATE)
	curl -L $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda

create-env: ## create conda environment
	if ${CONDA} env list | grep ${CONDA_ENV}; then \
	    mamba env update -n ${CONDA_ENV} -f environment.yml && \
	    echo "Activating new environment and installing R packages" && \
	    source ${ACTIVATE} ${CONDA_ENV} && \
	    R -e 'options(timeout=300); devtools::install_github("colomemaria/epiAneufinder", upgrade = "always"); q()'; \
	else \
	    ${CONDA} install -n base -y -c conda-forge mamba && \
	    source ${ACTIVATE} base && \
	    echo "Creating epianuefinder environment and installing dependencies" && \
	    mamba env create -f environment.yml && \
	    echo "Activating the environment and installing epianuefinder" && \
	    source ${ACTIVATE} ${CONDA_ENV} && \
	    R -e 'options(timeout=300); devtools::install_github("colomemaria/epiAneufinder", upgrade = "always"); q()'; \
	fi

.PHONY: create-env


help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
