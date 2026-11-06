#!/bin/bash

echo "run model..."
mpiexec -n 56 ./roms TWS.in > out.log 
echo "complete!"
