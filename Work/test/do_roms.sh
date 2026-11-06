#!/bin/bash

echo "run model..."
mpiexec -n 48 ./roms TWS2.in > out2.log 
echo "complete!"
