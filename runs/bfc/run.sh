docker run \
    -v ~/work/rcents/ec/runs/bfc:/work/bfc:ro \
    -v ~/work/data:/data:rw \
    -v ~/work/data/:/ncbi/public/:rw \
    rcents/im:bfc \
    scons --debug explain -C /data/ -f /work/bfc/rcebfc.py
