# "Disease-economy trade-offs under alternative pandemic control strategies" replication package

This repository contains data, code, and images for the paper "Disease-economy trade-offs under alternative pandemic control strategies". The model output was generated in R 3.6.3.

This repository is still being cleaned up and made into a turnkey experience; your patience is appreciated!


### How is this repository organized?

	* /Code contains all of the R code.

	* /Data contains the contact matrices.

	* /Figures contains all of the figures used in the paper.

### How do I recreate the figures and data files? (IN PROGRESS)

To run the code and generate the paper's figures and data:

	1. Make sure R and the necessary R packages are installed. You can find a list of necessary packages in the header of /Code/epicon_script.R

	2. Run "replicate_paper.sh". This will likely take some time. Note that the number of cores is set to 8 by default; you may want to adjust this to match your computing equipment and time requirements.

There are many parameters in "replicate_paper.sh" which are passed to subsequent scripts. These are described in the header of the file.

### What does the code in /bin do?

Broadly, there are three types of code files: R files which describe functions and algorithms (e.g., "functions.R"), R scripts which implement various algorithms and calculations, and shell scripts which call the R scripts.

### How is the data in /data organized?

#### Data inputs

#### Model outputs

