# docker run \
#     -v ~/work/rcents/ec/runs/bfc:/work/bfc:ro \
#     -v ~/work/data:/data:rw \
#     -v ~/work/data/:/ncbi/public/:rw \
#     rcents/im:bfc \
#     scons --debug explain -C /data/ -f /work/bfc/rcebfc.py

# docker run \
#     -v ~/work/rcents/ec/runs/bfc:/work/bfc:ro \
#     -v ~/work/data:/data:rw \
#     -v ~/work/data/:/ncbi/public/:rw \
#     rcents/im:bfc \
#     nextflow run -with-report /data/runs/rbfc.html -with-timeline /data/runs/tbfc.html -work-dir /data/runs/bfc /work/bfc/rbfc.nf

nextflow run rbfc.nf