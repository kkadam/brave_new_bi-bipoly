#!/bin/bash
### Set important parameters ###
run=m6
work_dir=/work/kkadam/scf_runs
greens=.true.

### Edit main.f90 to incorporate green's function array ###
sed -i -e '199d' main.f90
sed -i "199i\ have_green_funcs = $greens" main.f90


### Edit init file ??###
#sed -i -e '1,2d' init

make clean
make
#if ?$ then 


if [ -d $run ]; then
   rm -r $run
fi

mkdir $run
cp scf $run
cp init $run
cp runscf.h $run

sed -i -e '11d' bs
sed -i "11i\\$work_dir\/$run\/scf" bs

cp bs $run
if [ -d $work_dir\/$run ]; then
   echo "========================================================"
   echo "The run \"$run\" is ready in current working directory,"
   echo "because a directory with the same name is present in "
   echo "$work_dir "
   echo "========================================================"
else
   mv $run $work_dir
   echo "========================================================"
   echo "The run \"$run\" is ready in $work_dir!"
   echo "========================================================"
fi







