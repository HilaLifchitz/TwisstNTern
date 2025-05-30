> 💡 **Note for developers**  
> This is the draft for documentation and how command line implementation should look like (feel free to add/suggest more)
> All functions needs for command line implementation are provided in jupiter notebook, along this some testing.


# Running `twisstntern_simulate`

## Why simulate data?
The `twisstntern` module enables the implementation of a wide range of neutral demographic simulations. Simulating data allows researchers to explore how different demographic scenarios—such as population divergence, migration, and changes in population size—affect the distribution of genealogical topologies. These simulated patterns can be visualized using ternary plots, helping to build intuition about expected signals under neutrality and to compare against empirical genomic data.


## Usage

```
twisstntern_simulate -m <mode> -c <config_file> -o <output_dir>
```

### Required Arguments

* `-m <mode>`
  Simulation mode. Options:

  * `locus` – simulate a single locus
  * `chromosome` – simulate an entire chromosome

* `-c <config_file>`
  Path to the configuration file (see example below).

* `-o <output_dir>`
  Path to the directory where output files will be saved.

---

## Output

The following files will be generated in the specified output directory:

* `*.ts` — Tree sequence file containing the simulated genealogy.
* `ternary_plot.png` — Ternary plot visualizing ancestry contributions for your scenario.

---

## Example Configuration File

```
###########################################
# CONFIGURATION FILE FOR MSPRIME SIMULATION
###########################################

###########################################
STEP 1: Setting up demography
# Four-population demographic model
# NOTE: all simulations are haploid
# You can modify the values below or leave them as defaults

#=== SPLIT TIMES (in generations before present) ===#
# These define the timing of historical splits between populations.
# Larger values indicate older divergence events.

t1 = 1000         # Time of split between p1 and p2 (forms ancestor p12)
t2 = 5000         # Time of split between p12 and p3 (forms ancestor p123)
t3 = 10000        # Time of split between p123 and outgroup O

#=== EFFECTIVE POPULATION SIZES (diploid individuals) ===#
# Ne values for each extant and ancestral population.
# Affects genetic drift and coalescent times.

ne_p1 = 10000     # Population size of extant population 1
ne_p2 = 10000     # Population size of extant population 2
ne_p3 = 10000     # Population size of extant population 3
ne_p12 = 15000    # Size of ancestral population of p1 and p2
ne_p123 = 20000   # Size of ancestral population of p1, p2, and p3
ne_O = 30000      # Size of outgroup population

#=== MIGRATION RATES (per generation) ===#
# These rates are proportions: e.g. 0.001 means 0.1% of individuals replaced by migrants per generation.
# Define direction as "source > destination" (note that msprime defines in the reverse order)
# Default m = 0

# Migration between extant populations
m_p1>p2 = 0.001
m_p2>p1 = 0.001
m_p1>p3 = 0.0005
m_p3>p1 = 0.0005
m_p2>p3 = 0.0005
m_p3>p2 = 0.0005

# Migration between outgroup and extant populations
m_O>p1 = 0.0001
m_p1>O = 0.0001
m_O>p2 = 0.0001
m_p2>O = 0.0001
m_O>p3 = 0.0001
m_p3>O = 0.0001

# Migration involving ancestral populations
m_p12>p3 = 0.0002   # Migration from ancestor of p1/p2 to p3 (before t2)
m_p3>p12 = 0.0002
m_O>p123 = 0.00005  # Migration from outgroup to ancestral p123 (before t3)
m_p123>O = 0.00005

###########################################
STEP 2: Coalescent simulation
n_ind =  . Number of haploid individuals to include (default 20 per popuation)
[optional: mutation rate]
[optional: random seed]
[optional: give meaningful names to four populations]

#=== Simulation mode "locus" ===#
#msprime will simualte independent non-recombining loci with length 10k
n_loci =  . Number of loci / windows to simulate  (default 10000)

#=== Simulation mode "chromosome" ===#
msprime will simualte a chromosome of desired length and output tree sequences corresponding to breaks in recombination, genomic lenght of the loci simulated here will vary
rec_rate =    . Recombination rate per base per generation (e.g. 1e-8).

seq_length =  . Length of the chromosome to simulate (in base pairs).

###########################################

NOTE FOR DEVELOPERS: iteration through objects created by different modes will be different, will need to take care of it.

###########################################
# END OF CONFIGURATION
###########################################
```

