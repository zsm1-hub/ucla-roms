make -C Compile clean
rm -f Compile/roms roms

make -j1 \
  FFLAGS="-O3 -fp-model strict -traceback" \
  2>&1 | tee strict_make.log
