# r2r_init_bec2

This is a minimally modified UCLA `Tools-Roms/r2r_init.F`. It preserves the
original physical interpolation and adds standard BEC2 tracers found in the
parent restart file.

## Install and compile

Copy `r2r_init_bec2.F` into `Tools-Roms`, append
`r2r_init_bec2.Makefile.inc` to its `Makefile`, then run:

```bash
make clean
make FFLAGS="-O3 -qopenmp" LDFLAGS="-qopenmp" r2r_init_bec2
```

OpenMP must be present on both compile and link commands because shared
objects such as `roms_read_write.o` contain OpenMP symbols.

## Run

The arguments are identical to `r2r_init`:

```bash
./r2r_init_bec2 child_grid.nc theta_s theta_b hc N \
  parent_grid.nc parent_restart.nc rec
```

The output is `roms_init_bec2.nc`. The program always processes `temp` and
`salt`, then adds each standard BEC2 tracer that exists in the parent file.
Missing BEC2 variables are skipped, so the same executable works with BEC2
configurations that omit individual tracers.

The source currently recognizes the standard BEC2 spelling used by
the UCLA `bgc_real/code_check/tracers.opt`, with both `Alk` and `ALK` accepted.
The optional `NO2`, `N2`, and `N2O` tracers enabled by `Ncycle_SY` are also
recognized.
If the parent file uses different capitalization for other tracer names, edit
the `bec_name` block near the beginning of `r2r_init_bec2.F` to match `ncdump
-h` exactly; NetCDF variable names are case-sensitive.

## Verify

```bash
ncdump -h roms_init_bec2.nc | less
ncks -m roms_init_bec2.nc | egrep 'temp|salt|PO4|NO3|SiO3|NH4|Fe|O2|DIC|Alk|ALK|DOC|DON|DOFE|DOP|zooC|spC|diatC|diazC'
```

Check that every expected tracer has dimensions `(s_rho, eta_rho, xi_rho)`
as displayed by `ncdump`, and inspect surface and bottom maps after applying
the child `mask_rho`.
