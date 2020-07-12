#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************
#  This shell file gives examples of launching simulations for all 4 categories of workloads


#set this variable to NumCores in your cluster machine for faster sims
num_parallel_jobs=8

###########  HOW TO RUN JOBS?  ################

# The following line will launch sims for all workloads when you run ./doit.sh (comment it if you dont want it to) 

# ./runall.pl -s ../sim/predictor -w all -f  $num_parallel_jobs -d ../results/SEZNEC2014.08KB

./runall.pl -s ../sim/predictor -w random -f  $num_parallel_jobs -d ../results/results


###########  HOW TO GET STATS?  ################

# This scripts creates stats, after all the earlier jobs finish

./getdata.pl -w random -s MISPRED_PER_1K_INST -d ../results/results
./getdata.pl -w random -s BR_PRED_ACCURACY -d ../results/results

#./getdata.pl -w all -d ../results/SEZNEC2014.08KB

# To compare MPKI numbers against GSHARE for the provided benchmarks , uncomment this line 
# ./getdata.pl -w all -d ../results/MYRESULTS ../results/GSHARE.04KB  ../results/GSHARE.08KB ../results/GSHARE.16KB ../results/GSHARE.32KB



###########  PAST RUNS FOR 4KB, 8KB, and 16KB ###########

# The results are already in "../results" directory

#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.04KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.08KB"
#./runall.pl -f $num_parallel_jobs -d "../results/GSHARE.16KB"

################## GOOD LUCK! ##################
