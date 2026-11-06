#!/bin/bash

echo "run model..."
mpirun -n 8 ./roms USWC_wec_ana.in > out.log
echo "complete!"
