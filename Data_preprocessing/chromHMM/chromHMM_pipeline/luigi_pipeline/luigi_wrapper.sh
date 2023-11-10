export PYTHONPATH="/project2/lbarreiro/programs/bioluigi/:/project2/lbarreiro/programs/chromHMM/"
module load java

/project2/lbarreiro/programs/luigi/bin/luigi --module luigi_pipeline.chromHMM_l --local-scheduler $@
