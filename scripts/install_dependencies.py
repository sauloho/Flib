#!/usr/bin/env python

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""This script can be used to install the following Flib dependencies:
    * Blast v2.2.27
    * HH-suite v3.0
    * PSIPRED v4.0

It will also download the following databases:
    * HHblits (LATEST) [~39 Gb]
    * Blast (pdbaa)
    * Protein Data Bank (LATEST) [~94 Gb]

This script was written by Felix Simkovic, 2016.

===================================================
ANY USE OF THIS SCRIPT IS ENTIRELY AT YOUR OWN RISK
===================================================
"""

__author__ = "Felix Simkovic"
__date__ = "08 Aug 2016"
__version__ = 0.2

import argparse
import glob
import gzip
import logging
import os
import pip
import shutil
import sys
import tarfile
import tempfile

# We need wget all the way through and this is a good alternative 
# to normal command line invokes
try:
    import wget
except ImportError:
    pip.main(['install', 'wget'])
    import wget

# Setup some basic logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Not yet available for Windows so ignore
if sys.platform == "win32":
    logging.critical("Windows is not yet supported")

# ===============================================================================
# Show some command line options
# ===============================================================================

__PROGS = ['blast', 'hhsuite', 'psipred']
__DBS = ['blast', 'hhblits', 'pdb']

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-b', '--blast_db', type=str, default='pdbaa',
                    choices=['pdbaa'], 
                    help='Blast database [default: pdbaa]')
parser.add_argument('-t', '--hhblits_db', type=str, default='pdb70_04Aug16', 
                    choices=['pdb70_04Aug16', 'pfamA_29.0', 'uniprot20_2016_02'], 
                    help='HHblits database [default: pdb70_04Aug16]')
parser.add_argument('-f', '--flib_dir', type=str, default=os.getcwd(), help='Flib root directory [default: .]')
parser.add_argument('-i', '--install_dir', type=str, default=os.getcwd(), help='directory to install dependencies to [default: .]')
parser.add_argument('-o', '--overwrite', action="store_true", help='overwrite the current installations')
parser.add_argument('-p', '--progs', choices=__PROGS, default=__PROGS, nargs='+', help='programs to install [default = all]')
parser.add_argument('-d', '--databs', choices=__DBS, default=__DBS, nargs='+', help='databases to install [default = all]')
args = parser.parse_args()


# ===============================================================================
# Variable setting and processing of some options
# ===============================================================================

_BLAST_DB_NAME = args.blast_db
_HHBLITS_DB_NAME = args.hhblits_db
_DEPENDENCY_DIR = os.path.abspath(args.install_dir)
_DATABASE_DIR = os.path.join(_DEPENDENCY_DIR, 'databases')
_FLIB_DIR = os.path.abspath(args.flib_dir)
_OVERWRITE = True if args.overwrite else False
_PROGS_TO_INSTALL = args.progs
_DBS_TO_INSTALL = args.databs

# Make the installation directory if it doesn't exist
if not os.path.isdir(_DEPENDENCY_DIR): os.makedirs(_DEPENDENCY_DIR)

logging.info("""Installation options:

    * Programs to install: {0}
    * Databases to download: {1}

    * Overwrite: {2}

    * Flib directory: {3}
    * Database directory: {4}
    * Dependency directory: {5}
""".format(", ".join(_PROGS_TO_INSTALL), ", ".join(_DBS_TO_INSTALL), _OVERWRITE,
           _FLIB_DIR, _DATABASE_DIR, _DEPENDENCY_DIR)
)

# Some folder processing
if not os.path.isdir(_DATABASE_DIR):
    os.makedirs(_DATABASE_DIR)
if not os.path.isdir(_DEPENDENCY_DIR):
    os.makedirs(_DEPENDENCY_DIR)

# ===============================================================================
# Some functions to make the installation process easier
# ===============================================================================

def is_installed(folder, check_file):
    if os.path.isdir(folder) and check_file and os.path.isfile(check_file):
        return True
    else:
        return False

def do_install(folder, check_file=None, overwrite=False):
    """Simple check for installation of program based on folder"""
    if os.path.isdir(folder) and overwrite:
        return True
    elif is_installed(folder, check_file):
        return False
    else:
        return True

def print_overwrite(package):
    msg = "To re-install {0}, use \'--overwrite True\'".format(package)
    logging.info(msg)

def print_skipping(package):
    msg = "Not installing {0}".format(package)
    logging.info(msg)

# ===============================================================================
# Download and formatting of PDB database
# ===============================================================================

pdb_dir = os.path.join(_DATABASE_DIR, "pdb")
log_fname = tempfile.NamedTemporaryFile(delete=False).name
install_fname = os.path.join(pdb_dir, "install.ok")

if 'pdb' in _DBS_TO_INSTALL and do_install(pdb_dir, check_file=install_fname):
    logging.info("Creating a copy of the RCSB Protein Data Bank")

    if os.path.isdir(pdb_dir): shutil.rmtree(pdb_dir)

    rsync_cline = "rsync -rlpt -v -z --delete --port=33444 " \
                  "rsync.wwpdb.org::ftp/data/structures/divided/pdb/ {mirror_dir} " \
                  "> {log_dir} 2>&1"
    rsync_cline = rsync_cline.format(mirror_dir=pdb_dir, log_dir=log_fname)
    logging.info("Starting download - this might take some time ...")
    os.system(rsync_cline)

    # TODO: Determine differences and only extract new structures
    logging.info("Extracting PDB files to Flib readable tree")
    for folder in glob.glob(os.path.join(pdb_dir, "*")):
        for f in glob.glob(os.path.join(folder, "*")):
            name = os.path.basename(f).rsplit('.', 2)[0][3:]
            with gzip.open(f, 'rb') as f_in, open(os.path.join(pdb_dir, name+".pdb"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        shutil.rmtree(folder)

    os.unlink(log_fname)
    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Installation of RCSB Protein Data Bank completed")
    
elif 'pdb' in _DBS_TO_INSTALL:
    print_overwrite('PDB database')

elif is_installed(pdb_dir, install_fname):
    pass

else:
    print_skipping('PDB database')
    pdb_dir = ""

# ===============================================================================
# Download and extract the HHblits database
# ===============================================================================

_hhblits_db_dir = os.path.join(_DATABASE_DIR, "hhsuitedb", _HHBLITS_DB_NAME)
install_fname = os.path.join(_hhblits_db_dir, "install.ok")
hhblits_db = os.path.join(_hhblits_db_dir, _HHBLITS_DB_NAME)

if 'hhblits' in _DBS_TO_INSTALL and do_install(_hhblits_db_dir, check_file=install_fname):
    logging.info("Downloading the HHblits {0} database - please be very patient".format(_HHBLITS_DB_NAME))

    if os.path.isdir(_hhblits_db_dir): shutil.rmtree(_hhblits_db_dir)

    os.makedirs(_hhblits_db_dir)
    os.chdir(_hhblits_db_dir)
    tarball = wget.download("http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/{0}.tgz".format(_HHBLITS_DB_NAME))
    with tarfile.open(tarball) as tar: tar.extractall()
    if os.path.isfile(tarball) and os.path.isfile(hhblits_db): os.remove(tarball)
    os.chdir(_DEPENDENCY_DIR)
        
    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Download of HHblits {0} database completed".format(_HHBLITS_DB_NAME))

elif 'hhblits' in _DBS_TO_INSTALL:
    print_overwrite('HHblits database')

elif is_installed(_hhblits_db_dir, install_fname):
    pass

else:
    print_skipping('HHblits database')
    hhblits_db = ""

# ===============================================================================
# Download and extract the BLAST database
# ===============================================================================

_blast_db_dir = os.path.join(_DATABASE_DIR, "blastdb", _BLAST_DB_NAME)
install_fname = os.path.join(_blast_db_dir, "install.ok")
blast_db = os.path.join(_blast_db_dir, _BLAST_DB_NAME)

if 'blast' in _DBS_TO_INSTALL and do_install(_blast_db_dir, check_file=install_fname):
    logging.info("Downloading the Blast {0} database".format(_BLAST_DB_NAME))
        
    if os.path.isdir(_blast_db_dir): 
        shutil.rmtree(_blast_db_dir)

    os.makedirs(_blast_db_dir)
    os.chdir(_blast_db_dir)
    tarball = wget.download("ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz".format(_BLAST_DB_NAME))
    with tarfile.open(tarball) as tar: tar.extractall()
    if os.path.isfile(tarball): os.remove(tarball)
    os.chdir(_DEPENDENCY_DIR)

    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Download of Blast {0} database completed".format(_BLAST_DB_NAME))

elif 'blast' in _DBS_TO_INSTALL:
    print_overwrite('BLAST database')
 
elif is_installed(_blast_db_dir, install_fname):
    pass

else:
    print_skipping('BLAST database')
    blast_db = ""
    
# ===============================================================================
# Download and install the BLAST
# ===============================================================================

blast_dir = os.path.join(_DEPENDENCY_DIR, "blast")
install_fname = os.path.join(blast_dir, "install.ok")

if 'blast' in _PROGS_TO_INSTALL and do_install(blast_dir, check_file=install_fname):
    logging.info("Installing Blast legacy build") 
    
    if os.path.isdir(blast_dir): shutil.rmtree(blast_dir)

    if sys.platform == "linux" or sys.platform == "linux2":
        tarball = wget.download("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/blast-2.2.26-x64-linux.tar.gz")
    elif sys.platform == "darwin":
        tarball = wget.download("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/blast-2.2.26-universal-macosx.tar.gz")
    with tarfile.open(tarball) as tar: tar.extractall()
    if os.path.isfile(tarball): os.remove(tarball)
    shutil.move("blast-2.2.26", blast_dir)
    
    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Installation of Blast completed")

elif 'blast' in _PROGS_TO_INSTALL:
    print_overwrite('BLAST')

elif is_installed(blast_dir, install_fname):
    pass

else:
    print_skipping('BLAST')
    blast_dir = ""

# ===============================================================================
# Download and install PSIPRED
# ===============================================================================

psipred_dir = os.path.join(_DEPENDENCY_DIR, "psipred")
install_fname = os.path.join(psipred_dir, "install.ok")

if 'psipred' in _PROGS_TO_INSTALL and do_install(psipred_dir, check_file=install_fname):
    logging.info("Installing PSIPRED") 

    if os.path.isdir(psipred_dir): shutil.rmtree(psipred_dir)

    # tarball = wget.download("http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred.4.0.tar.gz")
    tarball = wget.download("https://github.com/psipred/psipred/archive/v4.0.tar.gz")
    with tarfile.open(tarball) as t_in: t_in.extractall(psipred_dir)
    if os.path.isfile(tarball): os.remove(tarball)
    for f in glob.glob(os.path.join(psipred_dir, "*")):
        for _f in glob.glob(os.path.join(f, "*")):
            shutil.move(_f, psipred_dir)
        shutil.rmtree(f)

    logging.info("Compiling PSIPRED")
    os.chdir(os.path.join(psipred_dir, 'src'))
    os.system('make')  # Compile from source
    for prog in ['chkparse', 'psipass2', 'psipred', 'seq2mtx']: 
        os.remove(os.path.join('..', 'bin', prog))
        shutil.move(prog, os.path.join('..', 'bin')) 
    os.chdir(_DEPENDENCY_DIR)

    runpsipred_lines = open(os.path.join(psipred_dir, "runpsipred"), "r").readlines()
    with open(os.path.join(psipred_dir, "runpsipred"), "w") as f_out:
        for line in runpsipred_lines:
            line = line.strip("\n")
            if line.startswith("set dbname ="):
                line = "set dbname = {0}".format(blast_db) 
            elif line.startswith("set ncbidir ="):
                line = "set ncbidir = {0}".format(os.path.join(blast_dir, "bin")) 
            elif line.startswith("set execdir ="):
                line = "set execdir = {0}".format(os.path.join(psipred_dir, "bin"))
            elif line.startswith("set datadir ="):
                line = "set datadir = {0}".format(os.path.join(psipred_dir, "data"))
            f_out.write(line + os.linesep)

    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Installation of PSIPRED completed")

elif 'psipred' in _PROGS_TO_INSTALL:
    print_overwrite('PSIPRED')

elif is_installed(psipred_dir, install_fname):
    pass

else:
    print_skipping('PSIPRED')
    psipred_dir = ""

# ===============================================================================
# Download and install the HHsuite
# ===============================================================================

hhsuite_tmp = os.path.join(_DEPENDENCY_DIR, "_hhsuite")
hhsuite_dir = os.path.join(_DEPENDENCY_DIR, "hhsuite")
install_fname = os.path.join(hhsuite_dir, "install.ok")

if 'hhsuite' in _PROGS_TO_INSTALL and do_install(hhsuite_dir, check_file=install_fname):
    logging.info("Installing HHsuite") 
    
    if os.path.isdir(hhsuite_dir): shutil.rmtree(hhsuite_dir)

    repo = pip.git.Git(url="git+https://github.com/soedinglab/hh-suite")
    repo.obtain(hhsuite_tmp)
    repo.update_submodules(hhsuite_tmp)
    
    # Build HH-suite
    logging.info("Compiling HHsuite")
    _build_dir = os.path.join(hhsuite_tmp, "build")
    os.makedirs(_build_dir)
    os.chdir(_build_dir)
    os.system("cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G \"Unix Makefiles\" -DCMAKE_INSTALL_PREFIX={0} ..".format(hhsuite_dir))
    os.system("make")
    os.system("make install")
    os.chdir(_DEPENDENCY_DIR)
    
    shutil.rmtree(hhsuite_tmp)
    with open(install_fname, 'w'): os.utime(install_fname, None)
    logging.info("Installation of HHsuite completed")

elif 'hhsuite' in _PROGS_TO_INSTALL:
    print_overwrite('HHsuite')

elif is_installed(hhsuite_dir, install_fname):
    pass

else:
    print_skipping('HHsuite')
    hhsuite_dir = ""

# ===============================================================================
# Modify the pipeline.sh script
# ===============================================================================

flib_script = os.path.join(_FLIB_DIR, "runflibpipeline")
pipeline_lines = [line.strip("\n") for line in open(flib_script, "r").readlines()]

logging.info("Adding dependency paths to the Flib \'runflibpipeline\' script")
with open(flib_script, "w") as f_out:
    for line in pipeline_lines:

        if line.startswith("export PSIPRED"):
            line = "export PSIPRED={0}".format(psipred_dir)

        elif line.startswith("export BLASTDB"):
            line = "export BLASTDB={0}".format(blast_db)

        elif line.startswith("export HHSUITE"):
            line = "export HHSUITE={0}".format(hhsuite_dir) 

        elif line.startswith("export HHLIB"):
            line = "export HHLIB={0}".format(hhsuite_dir) 

        elif line.startswith("export HHBLITSDB"):
            line = "export HHBLITSDB={0}".format(hhblits_db) 

        elif line.startswith("export FLIB"):
            line = "export FLIB={0}".format(_FLIB_DIR)

        elif line.startswith("export PDB"):
            line = "export PDB={0}".format(pdb_dir)

        f_out.write(line + os.linesep)

