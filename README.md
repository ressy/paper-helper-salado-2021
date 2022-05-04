# Data Gathering and Analysis for Salado 2021

Some scripts to aggregate data and reproduce outputs from:

Salado I, Fernández-Gil A, Vilà C, Leonard JA (2021) Automated genotyping of
microsatellite loci from feces with high throughput sequences. PLOS ONE 16(10):
e0258906. <https://doi.org/10.1371/journal.pone.0258906>

    git submodule update --init
    conda env update --file environment.yml
    conda activate paper-helper-salado-2021
    ./chiimp/install_linux.sh
    export PATH="$CONDA_PREFIX/lib/R/library/chiimp/exec:$PATH"
