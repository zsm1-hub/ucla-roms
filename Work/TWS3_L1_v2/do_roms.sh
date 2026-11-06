#!/bin/bash

echo "run model..."

nohup mpiexec -n 49 ./roms TWS.in > out.log 2>&1 &


echo "ROMS started in background"
echo "PID: $!"

