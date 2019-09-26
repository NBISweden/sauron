Bootstrap: docker
From: continuumio/miniconda3

%files
    environment.yml /tmp/environment.yml

%post
    /opt/conda/bin/conda env create -f /tmp/environment.yml
    echo "source activate $(head -1 /tmp/environment.yml | cut -d' ' -f2)" > ~/.bashrc

%environment
    PATH=/opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

%runscript
  exec "$@"
