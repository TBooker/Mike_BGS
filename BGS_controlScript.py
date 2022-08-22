import sys, argparse, re, os, tskit, pyslim, msprime, glob
from multiprocessing import Pool
import numpy as np
import pandas as pd
## Control script for Mike's BGS project

def analyseTrees(file_name, effective_size):
    # Load the tree sequence
    orig_ts = tskit.load(file_name)

    # Sprinkle on top the neutral mutations
    ts = msprime.sim_mutations( orig_ts,
    # Choose a neutral mutation rate to acheive neutral diversity of 0.005
                                rate = 0.005/(4*effective_size),
                                model = msprime.SLiMMutationModel(type = 0),
                                keep = True)

    # Set up window variables representing the distance from the site of deleterious mutations
    windows = np.array([0, 5000, 10000,11000])

    # Calculate nucleotide diversity in those windows
    d = ts.diversity(windows=windows)

    # grab diversity for the two windows:
    print(d[:2])

    # extract metadata for the simulation:
    meta_data_raw = re.sub(".trees", "", file_name).split("_")

    results = {meta_data_raw[n]:meta_data_raw[n+1]  for n in range(0, len(meta_data_raw),2)}

    results["window_1"] =  d[0]
    results["window_2"] =  d[1]

    return(results)


def slim_wrapper_function(list_of_args):



    command = " ".join([ list_of_args[0] ,
                        "-d","N="+str(list_of_args[1]),
                        "-d","Nes="+str(list_of_args[2]),
                        "-d","r="+str(list_of_args[3]),
                        "-d","rep="+str(list_of_args[5]),
                        list_of_args[4]])
#    print(command)
    file_name = 'Nes_'+ str(list_of_args[2]) + '_rep_'+ str(list_of_args[5])+ '_N_'+ str(list_of_args[1])+ '_r_'+ str(list_of_args[3])+ '.trees'
    os.system(command)
    print(file_name)
    return(file_name)

def main():

## Define command line args
    parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

    parser.add_argument("--population_size", "-N",
            required = True,
            dest = "pop_size",
            type = int,
            help = "The effective population size you want to simulate")

    parser.add_argument("--path_to_config","-p",
            required = True,
            dest = "path_to",
            type = str,
            help = "The location and name of the SLiM config script")

    parser.add_argument("--path_to_slim","-s",
            required = True,
            dest = "slim",
            type = str,
            help = "The full path to the SLiM build you want to use (assumes SLiM 4.x)")

    parser.add_argument("--recombination_rate","-r",
            required = True,
            dest = "rec",
            type = float,
            help = "The desired recombination rate between the functional sites and the focal neutral sites")

    parser.add_argument("--replicates",
            required = False,
            dest = "reps",
            type = int,
            help = "The number of replicates you want to run for each set of parameters",
            default = 50)

    parser.add_argument("--procs","-t",
            required = False,
            dest = "procs",
            type = int,
            help = "The number of threads to use",
            default = 1)

    parser.add_argument("--output","-o",
            required = True,
            dest = "output",
            type = str,
            help = "The name you want to give to the file")

    parser.add_argument("--simulate_only",
            dest = "simulate_only",
            action = "store_true",
            help = "This flag will skip the analysis step")

    parser.add_argument("--analysis_only",
            dest = "analysis_only",
            action = "store_true",
            help = "This flag will skip the simulate step and analyse all files matching *.trees in the current directory")


# Process command line args
    args = parser.parse_args()

    if args.analysis_only:
        output_trees_list = [g for g in glob.glob("*trees")]
    else:
    # The set of Ns values that you are going to analyse
        Ns_raw = np.array( [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 8, 10, 12, 15, 18, 20, 22, 24, 25, 28, 30, 50, 100, 200] )

    # convert Ns values to s for use in SLiM (add in the neutral case as well)
        selection_coefficients = np.concatenate([np.array([0]),
                                        Ns_raw/args.pop_size])
    # make a vector of replicates
        reps_vector = [args.reps*5] + [args.reps]*len(Ns_raw)

    # save a list of arguments for the map function
        argument_list = []

        for i in range(len(reps_vector)):
            tmp_list = [args.slim,
                        args.pop_size,
                        round(selection_coefficients[i]*args.pop_size, 5),
                        args.rec,
                        args.path_to ]
            for rep in range(reps_vector[i]):
                argument_list.append( tmp_list+[rep] )

    # Map the argument list to the slim wrapper function
        output_trees_list = []

        with Pool(args.procs) as p:
            tmp = p.map( slim_wrapper_function, argument_list)
            output_trees_list.append(tmp)

    if args.simulate_only:
        pass
    else:
    # Analyse each simulation and store the results
        results = [ analyseTrees(tmp_file, args.pop_size) for tmp_file in output_trees_list ]

    # Combine the results into a DF
        results_df = pd.DataFrame(results)

    # Save the final DF
        results_df.to_csv( args.output, index = False)


if __name__ == "__main__":
    main()
