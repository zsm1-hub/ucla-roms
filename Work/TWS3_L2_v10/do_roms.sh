#!/bin/bash

echo "run model..."

nohup mpiexec -n 50 ./roms TWS2.in > out.log 2>&1 &


echo "ROMS started in background"
echo "PID: $!"

