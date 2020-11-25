#!/bin/bash

#===================================================================#
# Description: Execute CFD-DEM parallel run
# Author: Wang Ding
#===================================================================#

# source CFDEM env vars
. ~/.bashrc
. $CFDEM_bashrc

# include functions
source $CFDEM_SRC_DIR/cfdem/etc/functions.sh

# define variables for mpirun
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
solverName="cfdemSolverIB-1.1.0"
numberOfProcs="4"
machineFileName="none" # yourMachinefileName
logPath=$casePath
logfileName="log_mpirun_$numberOfProcs_$solverName"

echo "Current case path:$casePath"

# block mesh
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
  echo "blockMesh: mesh was built before and using old mesh"
else
  echo "blockMesh: mesh needs to be built"
  cd $casePath/CFD
  blockMesh
fi

# call function to run a parallel CFD-DEM case
echo "Start mpirun -$numberOfProcs $solverName"
cfdemParallelRun $casePath $solverName $numberOfProcs $machineFileName $logPath $logfileName

# define variables for post-process
runOctave="false"
postProc="false"
liggghtsDumpFileName="dump.liggghts_run"

if [ $runOctave == "true" ]
then
  cd $casePath/CFD/octave
  octave --no-gui postproc.m
fi

if [ $postProc == "true" ]
then
  # get VTK data from liggghts dump file
  cd $casePath/DEM/post
  python -i $CFDEM_LPP_DIR/lpp.py $liggghtsDumpFileName
fi

#- keep terminal open (if started in new terminal)
echo "press Ctr+C kill process"
read
