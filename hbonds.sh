#!/bin/bash

bunzip2 ../DIP*
bunzip2 ../POS*

ln -s ../POSITION_CMD

# create list of all hydrogen bonding water species
./pos_collect_water_molecules input.collect

# create 1 and 2 donor hydrogen lists
mkdir sort_indexes
cd sort_indexes/
ln -s ../WATER.MOLINDEX
cp ../sort_MOLINDEX.sh .
./sort_MOLINDEX.sh
cd ..

# create 1 donor hydrogen DIP files
mkdir 01DH 2DH
cd 01DH/
ln -s ../sort_indexes/01DH_SORT_WATER.MOLINDEX
ln -s ../../DIPMOL_CMD
ln -s ../../DIPIND_CMD
cp ../input.buildmol_ind .
sed -i -e 's/FILE/01DH_SORT_WATER.MOLINDEX/g' input.buildmol_ind
ln -s ../../POSITION_CMD
../buildmol_ind input.buildmol_ind
rm DIPIND_CMD DIPMOL_CMD
ln -s DIPMOL_CMDtest DIPMOL_CMD
ln -s DIPIND_CMDtest DIPIND_CMD
cd ..

# create 2 donor hydrogens list
cd 2DH/
ln -s ../sort_indexes/2DH_SORT_WATER.MOLINDEX
ln -s ../../DIPMOL_CMD
ln -s ../../DIPIND_CMD
cp ../input.buildmol_ind .
sed -i -e 's/FILE/2DH_SORT_WATER.MOLINDEX/g' input.buildmol_ind
ln -s ../../POSITION_CMD
../buildmol_ind input.buildmol_ind
rm DIPIND_CMD DIPMOL_CMD
ln -s DIPMOL_CMDtest DIPMOL_CMD
ln -s DIPIND_CMDtest DIPIND_CMD
cd ..

