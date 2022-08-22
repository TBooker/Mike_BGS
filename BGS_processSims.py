import tskit, pyslim, msprime
import sys, re
import numpy as np

# Load the tree sequence
orig_ts = tskit.load(sys.argv[1])

# Sprinkle on top the neutral mutations
ts = msprime.sim_mutations( orig_ts,
# Choose a neutral mutation rate to acheive neutral diversity of 0.005
                            rate = 1.25e-7,
                            model = msprime.SLiMMutationModel(type = 0),
                            keep = True)

# Set up window variables representing the distance from the site of deleterious mutations
windows = np.array([0, 5000, 10000,11000])

# Calculate nucleotide diversity in those windows
d = ts.diversity(windows=windows)

# grab diversity for the two windows:
print(d[:2])

# extract metadata for the simulation:
print(sys.argv[1])

meta_data_raw = re.sub(".trees", "", sys.argv[1]).split("_")

results = {meta_data_raw[n]:meta_data_raw[n+1]  for n in range(0, len(meta_data_raw),2)}

results["window_1"] =  d[0]
results["window_2"] =  d[1]

print(results)
