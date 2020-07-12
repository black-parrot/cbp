#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************
#
#  This shell file gives examples of launching simulations for all 4 categories of workloads


#set this variable to NumCores in your cluster machine for faster sims
num_parallel_jobs=16

#../submissions/2014Unlimited/cbpUnlimited/predictor.cc
#../submissions/2014Unlimited/cpbUnlimited/predictor.h
#../submissions/AndreSeznecLimited/cbp64KB/predictor.cc
#../submissions/AndreSeznecLimited/cbp64KB/predictor.h
#../submissions/AndreSeznecLimited/cbp8KB/predictor.cc
#../submissions/AndreSeznecLimited/cbp8KB/predictor.h
#../submissions/AndreSeznecUnlimited/cbpUnlimited/predictor.cc
#../submissions/AndreSeznecUnlimited/cbpUnlimited/predictor.h
#../submissions/DanielJimenez1/cbp64KB/predictor.cc
#../submissions/DanielJimenez1/cbp64KB/predictor.h
#../submissions/DanielJimenez1/cbp8KB/predictor.cc
#../submissions/DanielJimenez1/cbp8KB/predictor.h
#../submissions/DanielJimenez1/cbpUnlimited/predictor.cc
#../submissions/DanielJimenez1/cbpUnlimited/predictor.h
#../submissions/DanielJimenez2/cbp64KB/predictor.cc
#../submissions/DanielJimenez2/cbp64KB/predictor.h
#../submissions/DanielJimenez2/cbp8KB/predictor.cc
#../submissions/DanielJimenez2/cbp8KB/predictor.h
#../submissions/DanielJimenez2/cbpUnlimited/predictor.cc
#../submissions/DanielJimenez2/cbpUnlimited/predictor.h
#../submissions/StephenPruett/cbp64KB/predictor.cc
#../submissions/StephenPruett/cbp64KB/predictor.h
#../submissions/StephenPruett/cbp8KB/predictor.cc
#../submissions/StephenPruett/cbp8KB/predictor.h




###########  HOW TO RUN JOBS?  ################

# The following line will launch sims for all workloads when you run ./doit.sh (comment it if you dont want it to) 
##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/2014Limited/cbp32KB/predictor.cc .
cp ../submissions/2014Limited/cbp32KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/2014_cbp32KB
./getdata.pl -w temp -d  ../results/2014_cbp32KB > ../results/2014_cbp32KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/2014Unlimited/cbpUnlimited/predictor.cc .
cp ../submissions/2014Unlimited/cbpUnlimited/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/2014_cbpUnl
./getdata.pl -w temp -d  ../results/2014_cbpUnl > ../results/2014_cbpUnl/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/AndreSeznecLimited/cbp8KB/predictor.cc .
cp ../submissions/AndreSeznecLimited/cbp8KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/AndreSeznec_cbp8KB
./getdata.pl -w temp -d ../results/AndreSeznec_cbp8KB > ../results/AndreSeznec_cbp8KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/AndreSeznecLimited/cbp64KB/predictor.cc .
cp ../submissions/AndreSeznecLimited/cbp64KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/AndreSeznec_cbp64KB
./getdata.pl -w temp -d ../results/AndreSeznec_cbp64KB > ../results/AndreSeznec_cbp64KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/AndreSeznecUnlimited/cbpUnlimited/predictor.cc .
cp ../submissions/AndreSeznecUnlimited/cbpUnlimited/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/AndreSeznec_cbpUnl
./getdata.pl -w temp -d  ../results/AndreSeznec_cbpUnl > ../results/AndreSeznec_cbpUnl/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez1/cbp8KB/predictor.cc .
cp ../submissions/DanielJimenez1/cbp8KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez1_cbp8KB
./getdata.pl -w temp -d ../results/DanielJimenez1_cbp8KB > ../results/DanielJimenez1_cbp8KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez1/cbp64KB/predictor.cc .
cp ../submissions/DanielJimenez1/cbp64KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez1_cbp64KB
./getdata.pl -w temp -d ../results/DanielJimenez1_cbp64KB > ../results/DanielJimenez1_cbp64KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez1/cbpUnlimited/predictor.cc .
cp ../submissions/DanielJimenez1/cbpUnlimited/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez1_cbpUnl
./getdata.pl -w temp -d ../results/DanielJimenez1_cbpUnl > ../results/DanielJimenez1_cbpUnl/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez2/cbp8KB/predictor.cc .
cp ../submissions/DanielJimenez2/cbp8KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez2_cbp8KB
./getdata.pl -w temp -d ../results/DanielJimenez2_cbp8KB > ../results/DanielJimenez2_cbp8KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez2/cbp64KB/predictor.cc .
cp ../submissions/DanielJimenez2/cbp64KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez2_cbp64KB
./getdata.pl -w temp -d ../results/DanielJimenez2_cbp64KB > ../results/DanielJimenez2_cbp64KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/DanielJimenez2/cbpUnlimited/predictor.cc .
cp ../submissions/DanielJimenez2/cbpUnlimited/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/DanielJimenez2_cbpUnl
./getdata.pl -w temp -d ../results/DanielJimenez2_cbpUnl > ../results/DanielJimenez2_cbpUnl/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/StephenPruett/cbp8KB/predictor.cc .
cp ../submissions/StephenPruett/cbp8KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/StephenPruett_cbp8KB
./getdata.pl -w temp -d ../results/StephenPruett_cbp8KB > ../results/StephenPruett_cbp8KB/rollup

##############################
cd ../sim
make clean
rm predictor.cc
rm predictor.h
cp ../submissions/StephenPruett/cbp64KB/predictor.cc .
cp ../submissions/StephenPruett/cbp64KB/predictor.h .
make
cd ../scripts
time ./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/StephenPruett_cbp64KB
./getdata.pl -w temp -d ../results/StephenPruett_cbp64KB > ../results/StephenPruett_cbp64KB/rollup

##############################


#./runall.pl -s ../sim/predictor -w temp -f  $num_parallel_jobs -d ../results/MYRESULTS


###########  HOW TO GET STATS?  ################

# This scripts creates stats, after all the earlier jobs finish

#./getdata.pl -w temp -d ../results/MYRESULTS
#./getdata.pl -w temp -d ../results/AWS_RESULTS_SEZNEC_EVAL

# To compare MPKI numbers against GSHARE for the provided benchmarks , uncomment this line 
# ./getdata.pl -w temp -d ../results/MYRESULTS ../results/GSHARE.04KB  ../results/GSHARE.08KB ../results/GSHARE.16KB ../results/GSHARE.32KB



###########  PAST RUNS FOR 4KB, 8KB, and 16KB ###########

# The results are already in "../results" directory

#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.04KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.08KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.16KB"

################## GOOD LUCK! ##################
