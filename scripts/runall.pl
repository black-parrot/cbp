#!/usr/bin/perl -w
#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************

require ( "./bench_list.pl");

#####################################
######### DEFAULT VARIABLES  ########
#####################################

$trace_dir = "../traces/";
$filetype  = ".bt9.trace.gz";
$wsuite    = "all"; 
$sim_exe   = "../sim/predictor";
$dest_dir  = "../results/MYRESULTS";
$debug     = 1;
$firewidth = 8; 

#####################################
######### USAGE OPTIONS      ########
#####################################

sub usage(){

$USAGE = "Usage:  '$0 <-option> '";
 
print(STDERR "$USAGE\n");
print(STDERR "\t-h                    : help -- print this menu. \n");
print(STDERR "\t-d <dest_dir>         : name of the result directory. \n");
print(STDERR "\t-w <workload/suite>   : workload suite from bench_list \n");
print(STDERR "\t-s <sim_exe>          : simulator executable \n");
print(STDERR "\t-dbg                  : debug \n");
print(STDERR "\t-f <val>              : firewidth, num of parallel simjobs to launch \n");
print(STDERR "\n");

exit(1);
}

######################################
########## PARSE COMMAND LINE ########
######################################
 
while (@ARGV) {
    $option = shift;
 
    if ($option eq "-dbg") {
        $debug = 1;
    }elsif ($option eq "-w") {
        $wsuite = shift;
    }elsif ($option eq "-s") {
        $sim_exe = shift;
    }elsif ($option eq "-d") {
        $dest_dir = shift;
    }elsif ($option eq "-f") {
        $firewidth = shift;
    }else{
	usage();
        die "Incorrect option ... Quitting\n";
    }
}

##########################################################             
# create full path to destination directory and copy sim
#$mysim = "$dest_dir"."/"."sim.bin";
$mysim =   $sim_exe;
system ("mkdir $dest_dir");
system ("chmod +x $mysim");
##########################################################
# get the suite names, and num_benchmarks from bench_list


die "No benchmark set '$wsuite' defined in bench_list.pl\n"
        unless $SUITES{$wsuite};

@workload_list = split(/\s+/, $SUITES{$wsuite});

$num_w = scalar @workload_list;


##########################################################

unless($debug){
    mkdir "$dest_dir";
    system ("ls $dest_dir");
#    system ("cp $sim_exe $mysim");
    system ("chmod +x $mysim");
}

##########################################################

#for($ii=0; $ii< $num_w; $ii++){
#    $bmkname   = $workload_list[$ii];
#    $outfile = $dest_dir. "/" . $bmkname . ".res";
#    $infile  = $trace_dir.$bmkname.$filetype;
#    $exe = "$mysim $infile > $outfile ";
#
#    #background all jobs except the last one (acts as barrier)
#    $exe .= " &" unless($ii%$firewidth==($firewidth-1));
#    print " Running $exe\n";
#    system("$exe") unless($debug);
#}
#

# Lets parallize the workload
$MAX_PROCESSES = $firewidth; #60;
my $numChildren=0; # number of forked children currently active

for($ii=0; $ii< $num_w; $ii++){
    $bmkname   = $workload_list[$ii];
    $outfile = $dest_dir. "/" . $bmkname . ".res";
    $infile  = $trace_dir.$bmkname.$filetype;
    $exe = "$mysim $infile > $outfile ";

    while ($numChildren >= $MAX_PROCESSES) {  # limit forked children to <= $processes
      my $deadKid = wait();
      $numChildren--;
    }
    $numChildren++;   # We're about to fork another kid
    my $pid = fork() and next; # Fork; the parent loops to do the next host and the child does the following:
    print "Running predictor on this trace $exe\n";
    system("$exe");
     exit;
}

  # Reap remaining kids
  while ($numChildren > 0) {
    my $deadKid = wait();
    $numChildren--;
  }

