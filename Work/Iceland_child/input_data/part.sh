echo "cp files here"
cp /leader/user/zsm/input-ucla/Iceland_child_era.nc ./
cp /leader/user/zsm/input-ucla/Iceland_child_bry.nc ./
cp /leader/user/zsm/input-ucla/Iceland_chd_initial_new.nc ./
mv Iceland_chd_initial_new.nc Iceland_chd_initial.nc 
cp /leader/user/zsm/input-ucla/Iceland_chgrd.nc ./

cp /leader/user/zsm/input-ucla/Iceland_child_tide.nc ./
echo "done"

./partit 6 6 Iceland_chgrd.nc  
./partit 6 6 Iceland_chd_initial.nc 
#./partit 6 6 Iceland_wnd.nc
#./partit 6 6 Iceland_tra.nc
#./partit 6 6 Iceland_rad.nc
#./partit 6 6 Iceland_prec.nc
./partit 6 6 Iceland_child_era.nc
./partit 6 6 Iceland_child_bry.nc 
./partit 6 6 Iceland_child_tide.nc
#./partit 6 6 Iceland_river.nc

