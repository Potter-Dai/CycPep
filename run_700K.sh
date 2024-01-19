#!/bin/bash

###Full Minimize###
mkdir min
cp APGfP.rsff2c.top inipep_frame.rst7 min/
cd min

/bin/cat > min.in<<mEOF
initial minimisation whole system
 &cntrl
  imin   = 1,
  irest  = 0,
  maxcyc = 40000,
  ncyc   = 20000,
  ntb    = 1,             
  ntr    = 0,            
  cut    = 9.0
  ntpr = 4000, ntwr = 4000,
 /
mEOF
CUDA_VISIBLE_DEVICES=6 pmemd.cuda -O -i min.in -o min.out -p APGfP.rsff2c.top -c inipep_frame.rst7 -r min.nrst 

############################
cd ..


mkdir heating
cp APGfP.rsff2c.top min/min.nrst heating/
cd heating
/bin/cat > heat.in<<hEOF
ACD-KTP : initial minimisation solvent
 &cntrl
  imin   = 0, 
  ig     = -1,  
  irest  = 0,   
  ntx    = 1,   
  ntb    = 1,   
  cut    = 9.0,
  ntr    = 0,   
  ntc    = 2,   
  ntf    = 2,   
  ntxo   = 2,   
  tempi  = 0.0,
  temp0  = 300.0,
  ntt    = 3,
  gamma_ln = 2.0,    
  nstlim = 100000, dt = 0.002,                     
  ntpr = 2000, ntwr = 2000, 
 /
hEOF
CUDA_VISIBLE_DEVICES=6 pmemd.cuda -AllowSmallBox -O -i heat.in -o heat.out -p APGfP.rsff2c.top -c min.nrst -r heat.nrst

############################
cd ..

mkdir equil
cp APGfP.rsff2c.top heating/heat.nrst equil/
cd equil

/bin/cat > equil.in<<eEOF
BCD-KTP : initial minimisation whole system
 &cntrl
  imin   = 0,
  ig     = -1,
  irest  = 0,
  ntx    = 1,
  ntb    = 2,
  pres0  = 1.0,
  ntp    = 1,
  cut    = 9.0,
  ntr    = 0,
  ntc    = 2,
  ntf    = 2,
  ntxo   = 2,
  tempi  = 300.0,
  temp0  = 300.0,
  ntt    = 3,
  gamma_ln = 2.0,
  nstlim = 500000, dt = 0.002,
  ntpr = 20000, ntwx = 20000, ntwr = 20000,
 /
eEOF
CUDA_VISIBLE_DEVICES=6 pmemd.cuda -AllowSmallBox -O -i equil.in -o equil.out -p APGfP.rsff2c.top -c heat.nrst -r equil.nrst -x equil.nc
#########################################
cd ..


for i in 1 2 3 4 5
do
mkdir replica_$i
cp APGfP.rsff2c.top equil/equil.nrst replica_$i
cd replica_$i

/bin/cat > replication.in<<EOF
replication simulation
 &cntrl
  imin   = 0,
  ig     = -1,
  irest  = 0,
  ntx    = 1,
  ntb    = 1,
  cut    = 9.0,
  ntc    = 2,
  ntf    = 2,
  ntxo   = 2,
  tempi  = 300.0,
  temp0  = 700.0,
  ntt    = 3,
  gamma_ln = 2.0,
  ioutfm = 1,
  nstlim = 300000000, dt = 0.002,
  ntpr = 25000, ntwx = 25000, ntwr = 25000,
 /
EOF
cd ..
done
# CUDA_VISIBLE_DEVICES=6 pmemd.cuda -AllowSmallBox -O -i replication.in -o replication_$i.out -p APGfP.rsff2c.top -c equil.nrst -r replica_$i.nrst -x replica_$i.nc
