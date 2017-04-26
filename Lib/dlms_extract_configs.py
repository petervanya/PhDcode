#!/usr/bin/env python
"""Usage:
    dlms_extract_configs.py <nf> [--levcfg <levcfg> --shift --nafion]
                                 [--histpath <hp>]

Extact config frames from DL_MESO HISTORY files, 
convert them to xyz files and place into Dump/ directory.
Need to have CONTROL file in the same directory.

Arguments:
    <nf>             Number of config frames requested [default: 10]

Options:
    --shift          true if shift of the box required
    --levcfg <lc>    Configuration level [default: 0]
    --nafion         Bead ordering "ABCWEP"
    --histpath <hp>  Path of the DL_MESO 'history_config.exe' file

pv278@cam.ac.uk, 21/06/16
"""
import os, sys, time
import subprocess
from subprocess import Popen, PIPE
from docopt import docopt


def get_timesteps():
    """Extract number of frames stored in HISTORY binaries
    from CONTROL file."""
    if not os.path.isfile("CONTROL"):
        sys.exit("No CONTROL file present.")
    with open("CONTROL") as f:
        for line in f:
            line = line.rstrip()
            if "traj" in line:
                Nstart, Nevery = list(map(int, line.split()[1:3]))
            if "steps" in line and not "equilibration" in line:
                Nsteps = int(line.split()[1])
    return (Nsteps - Nstart) // Nevery + 1


def get_cores():
    if not (os.path.isfile("HISTORY") or os.path.isfile("HISTORY000000")):
        sys.exit("No HISTORY files in the current directory.")
    pipe = Popen("ls -l HISTORY* | wc -l", stdout=PIPE, shell=True)
    return int(pipe.communicate()[0])


args = docopt(__doc__)
Nf = int(args["<nf>"])
levcfg = int(args["--levcfg"])
shift = "--shift" if args["--shift"] else ""
nafion = "--nafion" if args["--nafion"] else ""

Nhf = get_timesteps()
Nc = get_cores()
print("Cores: %i | Total num. frames: %i | Requested frames: %i" % \
      (Nc, Nhf, Nf))
            
if not os.path.exists("Dump"):
    os.makedirs("Dump")
if not os.path.exists("Dump_vel"):
    os.makedirs("Dump_vel")
if not os.path.exists("Dump_for"):
    os.makedirs("Dump_for")

hist_path = "~/sw/bin" if not args["--histpath"] else args["--histpath"]

history_exe = os.path.expanduser(hist_path + "/history_config.exe")
config2xyz = os.path.dirname(__file__) + "/config2xyz.py"
if not os.path.isfile(history_exe):
    sys.exit("Cannot find history_config.exe.")
if not os.path.isfile(config2xyz):
    sys.exit("Cannot find config2xyz.py.")

ti = time.time()
for i in range(Nhf-Nf+1, Nhf+1):
    cmd = "%s %i %i %i" % (history_exe, Nc, levcfg, i)
    subprocess.call(cmd, shell=True)

    cmd = "%s CONFIG.out %s %s" % (config2xyz, shift, nafion)
    subprocess.call(cmd, shell=True)

    os.rename("CONFIG.xyz", "Dump/dump_%04i.xyz" % i)
    if os.path.isfile("CONFIG.vel"):
        os.rename("CONFIG.vel",  "Dump_vel/dump_%03i.vel" % i)
    if os.path.isfile("CONFIG.for"):
        os.rename("CONFIG.for",  "Dump_for/dump_%03i.for" % i)
tf = time.time()
print("Total time: %.2f s." % (tf - ti))


