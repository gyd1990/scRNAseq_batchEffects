#!/usr/bin/env python
#
# generate 3 simulations
# then perform data quality control, call differentially distributed
# genes, compare results and do data characterization
# will write subdir into pwd
# the template directory and the generator script must be in pwd as well
#
# set parameters below

# simulation name
simName = "strongDD_strongBatch"
templateDir = "template"
generatorFile = "R/generate_datasets.R"

# dataset parameters
NBtchs = 3 # number of batches per group
NGns = 10000 # number of genes
NClls = 50 # number of cells per batch
NClls_CV = 0.4 # coeff of var for cells per batch

# poisson-beta parameters (all log2 transformed)
K = 4.2 # lambda controlling poisson distro (maximum mean)
K_CV = 0.25 # within batch coeff of var for K
K_CV2 = 0.08 # between batch coeff of var for K
A = 1.9 # alpha controlling beta distro
A_CV = 0.35 # within batch coeff of var for alpha
A_CV2 = 0.07 # between batch coeff of var for alpha
B = 3.5 # beta controlling beta distro
B_CV=0.4 # within batch coeff of var for beta
B_CV2 = 0.09 # between batch coeff of var for beta

# differential distributions
# 1/4 change in mean and var (k changes)
# 1/4 fixed mean, shape changes (a, b change equally)
# 1/4 shape and mean changes (a changes)
# 1/4 shape and mean changes (b changes)
ddP = 0.1 # proportion of differentially distributed genes (all types total)
ddM = 3 # factor for how much a param will change for d.d.

# hill parameters
KM = 4.5 # K_m controlling proportion of zeros
N = 3.5 # n (exp) controlling proportion of zeros
PZero_SD = 0.10 # std by which proportion of zeros varies




# libraries
import os
import shutil
import errno
import csv
import subprocess


# set some variables
mainDir = os.getcwd()
templateDir = os.path.abspath(templateDir)
generatorFile = os.path.abspath(generatorFile)
sims = ["sim1", "sim2", "sim3"]
paramsHead = ["NBtchs", "NGns", "NClls", "NClls_CV", "K", "K_CV", "K_CV2",
              "A", "A_CV", "A_CV2", "B", "B_CV", "B_CV2", "ddP", "ddM", "KM",
              "N", "PZero_SD"]
paramsVals = [NBtchs, NGns, NClls, NClls_CV, K, K_CV, K_CV2, A, A_CV, A_CV2, B,
              B_CV, B_CV2, ddP, ddM, KM, N, PZero_SD]


# defs

# copy all Rmd files
def copyRmd(index):
    Rmds = []
    for file in os.listdir("Rmd"):
        if file.endswith(".Rmd"): Rmds.append(file)
    for Rmd in Rmds:
        shutil.copy(os.path.join("Rmd", Rmd), Rmd)
    shutil.copy(os.path.join("Rmd", index), "index.Rmd")

# write yml file
def yml(name):
    with open("_bookdown.yml", "w") as ouf:
        ouf.write('rmd_files: ["index.Rmd"]\n')
        ouf.write('output_dir: "%s"\n' % name)

# remove all Rmd and other files
def cleanUp():
    Rmds = []
    for file in os.listdir(os.getcwd()):
        if file.endswith(".Rmd"): Rmds.append(file)
    [os.remove(Rmd) for Rmd in Rmds]
    os.remove("_bookdown.yml")
    shutil.rmtree("_bookdown_files")

# make book
def make(name, index):
    print "...working on %s" % name
    copyRmd(index)
    yml(name)
    opts = 'split_by = "chapter", config = list(toc = list(collapse = "section"))'
    cmd = 'bookdown::render_book("", bookdown::gitbook(%s))' % opts
    #subprocess.call(["Rscript", "-e", cmd])
    with open(os.devnull, "w") as ouf:
        subprocess.call(["Rscript", "-e", cmd], stdout = ouf, stderr = ouf)
    cleanUp()




# create dir tree and generate datasets
print "\n\nCreating %s..." % simName
os.mkdir(simName)
os.chdir(simName)
with open("params.csv", "w") as ouf:
    writer = csv.writer(ouf)
    writer.writerow(paramsHead)
    writer.writerow(paramsVals)
for sim in sims:
    shutil.copytree(templateDir, sim)
subprocess.call([generatorFile])
os.chdir(mainDir)


# analyse data
os.chdir(simName)
simDir = os.getcwd()
for sim in sims:
    os.chdir(sim)
    print "\nSimulation %s..." % sim
    make("DataQuality", "index_quali.Rmd")
    make("Benchmarking", "index_bench.Rmd")
    make("Characterization", "index_char.Rmd")
    os.chdir(simDir)


# done
print "\nDone\n"
