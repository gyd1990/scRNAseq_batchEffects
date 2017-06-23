#!/usr/bin/env python
#
# generate simulations over a grid of parameters
# then perform data quality control, call differentially distributed
# genes, compare and save results
# will write subdirs into pwd for each parameter setup
# the template directory and the generator script must be in pwd as well
#
# set parameters below
# then run "./run_gridSearch.py"

# Parameter Space
par_n = [50, 100, 200, 300] # average number of cells / batch
par_b = [0.5, 0.7, 1.0, 1.5] # factor for var of gene distros between batches
par_g = [0.5, 0.7, 1.0, 1.5] # factor for var of gene distros within batches
par_c = [0.5, 0.7, 1, 1.3] # factor for cell heterogeneity within batches

# important directories
templateDir = "template"
generatorFile = "R/generate_datasets.R"




# libraries
import os
import string
import shutil
import errno
import csv
import subprocess
import multiprocessing
from joblib import Parallel, delayed


# some variables
nproc = 24
mainDir = os.getcwd()
templateDir = os.path.abspath(templateDir)
generatorFile = os.path.abspath(generatorFile)


# defs

# copy directory recursively with
# files and subdirectories
def copy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

# create subdirs sim1:3 w/ template tree
# generate 3 datasets
# then do QC, Calling, Comparison
# in all 3 simulations
def simulate(sceneDir, mainDir, templ, gen):
    print "Working on: %s" % sceneDir
    os.chdir(sceneDir)
    pwd = os.getcwd()
    for sim in ["sim1", "sim2", "sim3"]:
        copy(templ, sim)
    subprocess.call([gen])
    cmd = 'rmarkdown::render("%s")'
    for sim in ["sim1", "sim2", "sim3"]:
        os.chdir(sim)
        with open(os.devnull, "w") as ouf:
            subprocess.call(["Rscript", "-e", cmd % "quali.Rmd"],
                            stdout = ouf, stderr = ouf)
            subprocess.call(["Rscript", "-e", cmd % "call.Rmd"],
                            stdout = ouf, stderr = ouf)
            subprocess.call(["Rscript", "-e", cmd % "compare.Rmd"],
                            stdout = ouf, stderr = ouf)
            subprocess.call(["Rscript", "-e", cmd % "estimate.Rmd"],
                            stdout = ouf, stderr = ouf)
        os.chdir(pwd)
    os.chdir(mainDir)




# bla bla
print "\nStarting grid search over parameter space.\n"
print "Provided n: %s" % ", ".join(str(i) for i in par_n)
print "Provided b: %s" % ", ".join(str(i) for i in par_b)
print "Provided g: %s" % ", ".join(str(i) for i in par_g)
print "Provided c: %s" % ", ".join(str(i) for i in par_c)


# create directories
# and write parameter combination into table
print "\n\nCreating trees..."
newDirs = []
for n in par_n:
    print "\n\nn = %d" % n
    for b in par_b:
        print "\nb = %.1f" % b
        for g in par_g:
            print "g = %.1f" % g
            for c in par_c:
                print "...creating c = %.1f" % c
                newDir = "n%d_b%.1f_g%.1f_c%.1f" % (n, b, g, c)
                newDirs.append(newDir)
                os.mkdir(newDir)
                os.chdir(newDir)
                with open("params.csv", "w") as ouf:
                    writer = csv.writer(ouf)
                    writer.writerow(["n", "b", "g", "c"])
                    writer.writerow([n, b, g, c])
                os.chdir(mainDir)


# do QC, Calling, Comparision
print "\n\n\nSimulating scenarios...\n"
Parallel(n_jobs = nproc)(delayed(simulate)(
    newDir, mainDir, templateDir, generatorFile) for newDir in newDirs)


# done
print "\n\nDone.\n"
