'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 1.0
CREATED: 4-10-2021
REVISED:
DESCRIPTION: Automatizes the assembly protocol
INPUT (by order):
OUTPUT:
USAGE:
    Execute the script in the directory where the service will be performed.

REQUIREMENTS:
    -Python 3.6 or higher (f-strings)
    -os
    -requests
    -datetime
    -argparse

DISCLAIMER:
TO DO:
    -Fix the lablog format in the offline mode (seems to be tabbed fsr)
================================================================
END_OF_HEADER
================================================================
'''
# Imports

import os
import requests
import datetime
import argparse

# Functions
def download_lablog(url):
    # Downloads the lablog from the repository and writes it to a file

    request = requests.get(url)
    
    with open("lablog","w") as outfile:
        outfile.write(request.text)
    return

# Argument parser
parser = argparse.ArgumentParser(description="Prepare the lablogs for an assembly service in BU-ISCIII")
parser.add_argument("--offline-mode", 
                     default=False, 
                     action="store_true", 
                     dest="offline", 
                     help="Create the lablogs through offline hardcoding (might not be updated)")
args = parser.parse_args()

# Basic scaffolding of the service
os.mkdir("ANALYSIS")
os.mkdir("DOC")
os.mkdir("RAW")
os.mkdir("REFERENCE")
os.mkdir("RESULTS")
os.mkdir("TMP")

# Change the directory and download the first lablog
os.chdir("ANALYSIS")

if args.offline:
    with open("lablog","w") as outfile:
        outfile.write("cd 00-reads; cat ../samples_id.txt | xargs -I % echo \"ln -s ../../RAW/%_*R1*.fastq.gz %_R1.fastq.gz\" | bash; cd -\n\
                       cd 00-reads; cat ../samples_id.txt | xargs -I % echo \"ln -s ../../RAW/%_*R2*.fastq.gz %_R2.fastq.gz\" | bash; cd -"
                      )
else:
    url = "https://raw.githubusercontent.com/GuilleGorines/BU-ISCIII_services/main/ASSEMBLY/lablog"
    download_lablog(url)


directory = "00-reads"
os.mkdir(directory)

# Get time for the directory name
date = datetime.date.today()
date = date.strftime("%Y%m%d")

analysis_dir = f"{date}_ANALYSIS01_ASSEMBLY"

os.mkdir(analysis_dir)
os.chdir(analysis_dir)

if args.offline:
    with open("lablog","w") as outfile:
        outfile.write(
            "ln -s ../00-reads .\n \
             ln -s ../samples_id.txt .\n \
             ln -s /processing_Data/bioinformatics/pipelines/bacterial_qc .\n \
             #conda activate nextflow \n \
             echo 'nextflow run /processing_Data/bioinformatics/pipelines/bacterial_assembly-nf/main.nf -bg -resume --reads \"00-reads/*_R{1,2}.fastq.gz\" --fasta ../../REFERENCES/GCF_002072775.2_ASM207277v2_genomic.fna --gtf ../../REFERENCES/GCF_002072775.2_ASM207277v2_genomic.gff --outdir 03-assembly -profile hpc_isciii' > _01_nf_assembly.sh \n  \
             #nohup bash _01_nf_assembly.sh &> $(date '+%Y%m%d')_assembly01.log &"
             )

else:
    url = "https://raw.githubusercontent.com/GuilleGorines/BU-ISCIII_services/main/ASSEMBLY/ANALYSIS_ASSEMBLY/lablog"
    download_lablog(url)

# get the parent dir to go back to
parent_dir = os.getcwd()

# Preprocessing
directory = "01-preprocessing"
os.mkdir(directory)
os.chdir(directory)

if args.offline:
    with open("lablog","w") as outfile:
        outfile.write(
            "cat ../samples_id.txt | while read in; do echo \"mkdir $in;qsub -V -b y -j y -cwd -N TRIMMOMATIC.$in -q all.q@obelix03 -q all.q@obelix09,all.q@obelix10,all.q@obelix11 -pe openmp 10 java -jar -Djava.io.tmpdir=../../TMP/ /opt/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 -phred33 ../00-reads/\"$in\"_R1.fastq.gz ../00-reads/\"$in\"_R2.fastq.gz $in/\"$in\"_R1_filtered.fastq $in/\"$in\"_R1_unpaired.fastq $in/\"$in\"_R2_filtered.fastq $in/\"$in\"_R2_unpaired.fastq ILLUMINACLIP:/opt/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50\"; done >> _01_preprocess.sh"
            )

else:
    url = "https://raw.githubusercontent.com/GuilleGorines/BU-ISCIII_services/main/ASSEMBLY/ANALYSIS_ASSEMBLY/01-preprocessing/lablog"
    download_lablog(url)
os.mkdir("logs")
os.chdir(parent_dir)

# Kmerfinder
directory = "02-kmerfinder"
os.mkdir(directory)
os.chdir(directory)

if args.offline:
    with open("lablog","w") as outfile:
        outfile.write(
            "#module load singularity/singularity-2.6.0\n\n \
            kmerFinder_DB=/processing_Data/bioinformatics/references/kmerfinder/20190108_stable_dirs/bacteria\n \
            kmerFinder_path=/processing_Data/bioinformatics/pipelines/kmerfinder_v3.0.simg \n\n \
            cat ../samples_id.txt | xargs -I % echo \"qsub -V -b y -j y -cwd -N KMERFINDER -q all.q@obelix05 -q all.q@obelix04 singularity run --bind $PWD:/media --bind $kmerFinder_DB:/mnt --bind ../01-preprocessing:/workdir $kmerFinder_path -i /workdir/%/%_R1_filtered.fastq.gz /workdir/%/%_R2_filtered.fastq.gz -o /media/% -db /mnt/bacteria.ATG -tax /mnt/bacteria.name -x\" > _01_kmerfinder.sh \n\n \
            echo \"cat ../samples_id.txt | xargs -I % awk '{FS=\\\"\\\t\\\"} NR==2 {print \$1}' %/results.txt | awk '{count[\$0]++} END{for (i in count) {print count[i], i}}' | sort -nr\" > _02_find_common_reference.sh \n \
            echo \"bash _02_find_common_reference.sh | head -n1 | tr ' ' '\t' | cut -f2 | while read in; do bash /processing_Data/bioinformatics/references/bacteria/download_reference.sh \${in} ../../../REFERENCES/; done\" > _03_download_reference.sh"
            )
else:   
    url = "https://raw.githubusercontent.com/GuilleGorines/BU-ISCIII_services/main/ASSEMBLY/ANALYSIS_ASSEMBLY/02-kmerfinder/lablog"
    download_lablog(url)

os.mkdir("logs")
os.chdir(parent_dir)

# Kmerfinder statistics
directory = "99-stats"
os.mkdir(directory)
os.chdir(directory)

if args.offline:
    with open("lablog","w") as outfile:
        outfile.write(
            "python3 ../bacterial_qc/parse_kmerfinder.py --path ../02-kmerfinder --output_bn kmerfinder.bn --output_csv kmerfinder.csv"
            )
else:
    url = "https://raw.githubusercontent.com/GuilleGorines/BU-ISCIII_services/main/ASSEMBLY/ANALYSIS_ASSEMBLY/99-stats/lablog"
    download_lablog(url)

print("Folders are now ready for the assembly service!")