#echo "cp files here"
#cp /leader/user/zsm/input-ucla/Iceland_* ./
#echo "done"

./partit 6 12 Iceland_roms_grid.nc 
./partit 6 12 Iceland_initial_field.nc
#./partit 6 6 Iceland_wnd.nc
#./partit 6 6 Iceland_tra.nc
#./partit 6 6 Iceland_rad.nc
#./partit 6 6 Iceland_prec.nc
./partit 6 12 Iceland_era.nc
./partit 6 12 Iceland_bounary_field_200001.nc
./partit 6 12 Iceland_bounary_field_200002.nc
./partit 6 12 Iceland_bounary_field_200003.nc
./partit 6 12 Iceland_bounary_field_200004.nc
./partit 6 12 Iceland_tidal_forcing.nc
./partit 6 12 Iceland_river.nc
