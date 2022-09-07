


# SSMPG 2022
Repository for [Software and Statistical Methods for Population Genetics (SSMPG 2022)](https://ssmpg.sciencesconf.org/) (Aussois, September 19-23 2022)


##  1. Install software

### Install R and Rstudio
To participate in the practical sessions, bring your own laptop and install [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/), an integrated development environment (IDE) for R.

### Install R packages (LEA, gradientForest, vegan)
To install R packages for the data analyses, copy and paste the following pieces of code in the R session

```r
#Install R packages for SSMPG 2022

#Install packages from github
install.packages("devtools")

#Package LEA (development/latest version) 
devtools::install_github("bcm-uga/LEA")

#Package lfmm (also in LEA as LEA::lfmm2)
install.packages("lfmm")

#Package gradientForest
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

#Package vegan for RDA
install.packages("vegan")

#Package q-value for controlling FDR
#Try https:// or http:// 
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")



```


### Install BAYPASS (to be updated)

Download the archive for the latest stable version (2.3) from http://www1.montpellier.inra.fr/CBGP/software/baypass/ or directly via the following command run on a terminal:
```
wget http://www1.montpellier.inra.fr/CBGP/software/baypass/files/baypass_2.3.tar.gz
```
Extract the archive, *e.g.*, from a terminal:
```
tar -zxvf baypass_2.3.tar.gz
```
The source files are to be found in the *src* subdirectory. BayPass is coded in Fortran90 and can therefore be compiled for any system supporting a Fortran90 compiler using the provided Makefile. This Makefile is designed to work with either the free compiler *gfortran* (if not already installed in your system, binaries are available at https://gcc.gnu.org/wiki/GFortranBinaries and are easy to install for most Windows, Mac and Linux OS versions) or the commercial *ifort* intel Fortran compiler that is now freely available (for non commercial use) as part of the *oneAPI* Toolkit (see installation instruction for Windows, MacOS and Linux system). 
BayPass also uses OpenMP (http://openmp.org/wp/) to implement multithreading, which allows parallel calculation on computer systems that have multiple CPUs or CPUs with multiple cores. Users thus have to make sure that the corresponding libraries are installed (which is usually the case, on Linux OS or following compiler installation previously described). The following instructions run within the *src* subdirectory allows to compile the code and to produce a binary:
* using the *gfortran* free compiler (the command should automatically produce an executable called *g_baypass*):
```
make clean all FC=gfortran
```
* using the *ifort* intel Fortran compiler (the command should automatically produce an executable called *i_baypass*):
```
make clean all FC=ifort 
```
> Note: Under Linux (or MacOS), before the first use, make sure to give appropriate execution rights to the program. For instance you may run:
>```chmod +x baypass```


##  2. Download simulated datasets

Attendees will have the opportunity to explore methods by studying simulated data, discussing best practices and methodological weaknesses of these studied techniques. Practical sessions will favor interactions among participants in a collaborative spirit. 

The first data set is for the **whooly marmot** populations. The whooly marmot is an emblematic species that lived in Thibaut's computer with GENEPOP DNA code and survived abrupt change in environmental conditions. Thibaut recorded the fitness loss of whooly marmot populations after this change and the goal is to be as close as possible to the ground-thruth measure of fitness loss with genetic offset measures.   The data consist of a matrix of genotype for $n = $ individuals genotyped at $L = $ loci. x environmental variables have been measured before and after the abrupt change in conditions. 

The second data set is for the **Osuah tree** populations. Osuah (pronounced Aussuah) in a small village in the French alps with virtual tree species that grow in computers only. The Osuah tree is an emblematic species of Clement's computer with SLIM DNA code that also survived an abrupt environmental change. Clement recorded the survival probability for each populations before and after change. The goal is to be as close as possible to Clement's ground-thruth measure of fitness loss with the genetic offset measures.  The data consist of a matrix of genotype for $n = 300$ individuals genotyped at $L = 2333$ loci. Four environmental variables have been measured before and after the abrupt change in conditions. 


## 3. Create research groups and submit your paper

During data analysis sessions,  participants are encouraged to form teams. A team can be composed of 1 to 5 participants. Teams will collectively report a synthesis of their analysis in public (last day), and send three files for each data analysis 

* a list of candidate loci detected by their prefered GEA method or combination of methods
* a list of offset values obtained from their prefered offset method (or combination of methods)
* a short README file that explains the choices made. 

Eeach team will be asked to present 2-3 slides for each challenge.

## 4. Evaluation

Thibaut and Clement will reveal the truth about their simulation and what could be inferred from genotypes and environmental data.  Don't worry. Everyone wins! 
