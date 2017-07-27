All models were tested using the iSAT tool which can be found, along with documentation, here:

http://www.avacs.org/tools/isatode/

We include a version of iSAT within this tool evaluation. In our experiments we use the version r2806.

The server our evaluation was tested upon runs Python 2.7.13 although we speculate most versions of Python 2.7 will be sufficient. 

The kinetics of the CRN generated from our program can be simulated using the free tool Microsoft Visual GEC. This can give the reader some insight into the validity of our results.

http://biology.azurewebsites.net/gec/


--------------------------
1. Introduction
--------------------------


We evaluate the correctness and performance of our parameter synthesis method on three case studies: bellshape, phosphorelay and poisson explained within the paper “Syntax-Guided Optimal Synthesis for Chemical Reaction Networks”. All experiments were run on a server with a Intel Xeon CPU E5645 @2.40GHz processor (using a single core) and 24GB @1333MHz RAM.

Each individual subfolder included contains a python file that can be used to call the SAT-solver in order to run the experiments and results obtained within the paper - specifically the experiments from the experimental evaluation. Also included is a python script which can be used for individual calls to the solver on each model.
 


--------------------------
2. Running the experiments described within the paper.
--------------------------

Step 1. In order to generate the outputs from the solver please use the python file included within each folder. This is individual for each experiments

3.


--------------------------
Bellshape
--------------------------

The cost constraints are labelled under “— cost constraint” within bellshape.hys and can be adjusted accordingly.

The interpretation of discrete parameters is also commented to guide on how to interpret output results.

An example call for the bellshape model would be:

./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i superpoisson10.hys --ode-opts=--continue-after-not-reaching-horizon --prabs=0.001 --msw=0.02 --max-depth 2


Where --prabs is the delta precision
	  --msw is the minimum splitting width (used to control rate precision) 0.01 can also be used for higher rate precision (seen in example candidate solution files). Change the value 50 in lines:
	INIT… 	  dK_dt = (Kd - K)*50;
	TRANS … 	(dK_dt = (Kd - K)*50);
to 200+ for increase in performance with lower msw.

	  --continue-after-not-reaching-horizon is used to force the solver to continue even if the muller bracketing method fails (it will then switch to another enclosure method)
	  --max-depth refers to the number of unwindings of iSAT.


--------------------------
Phosphorelay
--------------------------

An example call would be:

nohup ./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i phosphorelaywithcovariances.hys --ode-opts=--continue-after-not-reaching-horizon --prabs=0.00001 --msw=0.01  --max-depth 2 > phospoutput.txt

where nohup allows the program to continue even if disconnected from the server. This example would be for precision 10-5

or

nohup ./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i phosphorelaywithoutcovariances.hys --ode-opts=--continue-after-not-reaching-horizon --prabs=0.00001 --msw=0.01  --max-depth 2 > phospoutput.txt

without the covariances.

--------------------------
Super Poisson
--------------------------


	print "Welcome to the program designed to run the Experimental Evaluation for the superpoisson model described within the paper: Syntax-Guided Optimal Synthesis for Chemical Reaction Networks"
	print "Please enter a number and press enter for one of the following options:"
	print '1) Interval Experiment - A procedure that automates calls to the solver for each superpoisson model where each model increases the real intervals as described in the paper (table3).' 
	print '2) Custom Experiment - If you do not wish to perform all experiments above you can generate a custom experiment for a specific precision or model. '
	print '3) Construct ResultsSummary- Construct a runtime / satisfiability summary for every file within the results folder folder'
	print '4) Custom Result- Given a filename the CRN is constructed and assembled from the file with a summary of runtime and satisifiability.'

An example call would be:

./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i superpoisson.hys --ode-opts=--continue-after-not-reaching-horizon --prabs=0.001 --msw=0.01 --max-depth 2







