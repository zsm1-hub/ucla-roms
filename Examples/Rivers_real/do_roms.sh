#!/bin/bash

echo "run model..."
mpiexec -n 18 ./roms rivers.in
echo "complete!"
