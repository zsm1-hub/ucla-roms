#!/bin/bash

echo "run model..."
mpiexec -n 50 ./roms TWS.in > out.log 
echo "complete!"
