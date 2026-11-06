# Add these lines beside the existing R2R_INIT/r2r_init target.
# This target reuses every original r2r_init object except the main program.

R2R_INIT_BEC2 = r2r_init_bec2.o r2r_interp_init.o r2r_subs.o r2r_rotate.o \
                r2r_vert_interp.o set_depth.o etch_into_land.o spln2d.o \
                tiling.o roms_read_write.o read_string.o lenstr.o

r2r_init_bec2: $(R2R_INIT_BEC2)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o r2r_init_bec2 $(R2R_INIT_BEC2) $(LCDF)
