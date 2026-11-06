#!/bin/bash

echo "run model..."
mpiexec -n 36 ./roms Iceland.in    
echo "complete!"
