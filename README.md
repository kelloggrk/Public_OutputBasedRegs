# Public_OutputBasedRegs
## Repository for "Output and Attribute-Based Carbon Regulation Under Uncertainty", by Ryan Kellogg

### Organization
- All of the code for this project is stored in this repository. To run the code, you need to clone this repo to your local machine (or copy the files directly, preseving the subfolder structure).
  * The local folder holding the code files is referred to as `repodir` below.
- All raw data files are available [here](https://www.dropbox.com/sh/bbumsxyrngrcmqm/AACTPSZy0i1YwZ6g1ZjaSs3fa?dl=0) in the `rawdata` folder.
  * This link also includes intermediate (`intdata`) and output (`output`) data folders that are populated with the output from the Stata scripts.
  * The `rawdata` subfolders contain README.md files with additional info on how each raw data file was originally accessed.
  * Users should download the `rawdata`, `intdata`, and `output` folders together into a single folder. This folder is referred to as `dropbox` below.


### Stata
- Code has been verifed to run on Windows OS 64-bit Stata SE v15. No extra package installations are required.
- To run the Stata code, you need a file called `globals.do` stored in your local root OutputBasedRegs repo folder. (`globals.do` is .gitignored)
    - `globals.do` should look like the below, pointing to your own directories:
```
global repodir = "C:/Work/OutputBasedRegs"
global dropbox = "C:/Users/Ryan Kellogg/Dropbox/OutputBasedRegs"
```

### Matlab
- Code has been verified to run on Windows OS Matlab R2018b. No extra package installations are required.
- To run the Matlab code, you need a file called `globals.m` stored in your local root OutputBasedRegs repo folder. (`globals.m` is .gitignored)
    - `globals.m` should look like the below, pointing to your own directories:
```
repodir = 'C:/Work/OutputBasedRegs';
dropbox = 'C:/Users/Ryan Kellogg/Dropbox/OutputBasedRegs';
```
- mainscript.m is the front-end script that sets up the model parameters, instantiates the model objects, calls the methods that generate the results, and then generates the figures and welfare results. 
- The model objects are defined in indexmodel.m and indexmodelF.m. indexmodel.m is the superclass and only allows uncertainty in eta. indexmodelF.m is a subclass that allows uncertainty in eta and F.
- If you are unfamiliar with object-oriented programming in Matlab, see resources available [here](https://www.mathworks.com/discovery/object-oriented-programming.html)


### LaTex
- A full LaTex build is required to compile the paper.


### Batch script
- The file `makepaper.sh` will run all Stata, Matlab, and LaTex code in order, resulting in the final compiled paper. It must be called from a bash shell.
  * The root directory definitions in this script for the locally cloned repository and the dropbox data need to point to your own directories.
  * See the README header in `makepaper.sh` for guidelines on running the script.

