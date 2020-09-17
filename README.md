# YeastHybridModelingFramework

Integration of a Boolean nutrient signaling network with an enzyme-constrained model of *S. cerevisiae*'s metabolism.

This repository contains all the necessary data, model files and scripts for reproducing the results on the publication **"A novel yeast hybrid modeling framework integrating Boolean and enzyme-constrained networks enables exploration of the interplay between signaling and metabolism"** available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.09.11.290817v1).

- Abstract

The interplay between nutrient-induced signaling and metabolism plays an important role in maintaining homeostasis and malfunction of this interplay has been implicated in many different human diseases such as obesity, type 2 diabetes, cancer and neurological disorders. Therefore, unravelling the role of nutrients as signaling molecules and metabolites as well as their interconnectivity may provide a deeper understanding of how these conditions occur. Both signalling and metabolism have been extensively studied using various systems biology approaches. However, they are mainly studied individually and in addition current models lack both the complexity of the dynamics and the effects of the crosstalk in the signaling system. To gain a better understanding of the interconnectivity between nutrient signaling and metabolism in Eukaryotes, we developed a hybrid model by combining Boolean, describing the signalling layer, and enzyme constrained models  accounting for metabolism using a regulatory network as a link for the yeast Saccharomyces cerevisiae. The model was capable of reproducing the regulatory effects that are associated with the Crabtree effect and glucose repression. We show that using this methodology one can investigate intrinsically different systems, such as signaling and metabolism, in the same model and gain insight into how the interplay between them can have non-trivial effects by showing a connection between Snf1 signaling and chronological lifespan by the regulation of NDE and NDI usage in respiring conditions. In addition, the model showed that during fermentation, enzyme utilization is the more important factor governing the protein allocation, while in low glucose conditions robustness and control is prioritized. 

- Reference:  
>A novel yeast hybrid modeling framework integrating Boolean and enzyme-constrained networks enables exploration of the interplay between signaling and metabolism.
Linnea Österberg, Iván Domenzain, Julia Münch, Jens Nielsen, Stefan Hohmann, Marija Cvijovic
bioRxiv; doi: [10.1101/2020.09.11.290817](https://doi.org/10.1101/2020.09.11.290817)

- Last update: 2020-09-17

This repository is administered by [@linoste](https://github.com/linoste), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation

### Required software:
* A functional Matlab installation (MATLAB_2017b or higher).
* [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB.
* libSBML MATLAB API ([version 5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).
* [Gurobi Optimizer for MATLAB](http://www.gurobi.com/registration/download-reg).

### Installation Instructions
* Clone the [master](https://github.com/cvijoviclab/YeastHybridModelingFramework) branch from [cvijoviclab GitHub](https://github.com/cvijoviclab).
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).

## Development Guidelines

Anybody is welcome to contribute to the development of this modeling and simulation Toolbox, but please abide by the following guidelines.

Each function should start with a commented section describing the function and explaining the parameters. Existing functions can clarify what style should be used.
### Bugfixes, new features and functions
* For any development, whether bugfixes, introducing new functions or new/updated features for existing functions: make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`.
* Make commits to this branch while developing. Aim for backwards compatibility, and try to avoid very new MATLAB functions when possible, to accommodate users with older MATLAB versions.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries (model `MATLAB` structures), UniProt databases, physiology and protemics data files, etc.
* `doc`: updating documentation (in `doc` folder) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: exportModel additional export to YAML
chore: update UniProt database for CENPK113-7D
fix: optimizeProb parsing results from Gurobi
```
More detailed explanation or comments can be left in the commit description.

### Pull request
* No changes should be directly commited to the `master` or `devel` branches. Commits are made to side-branches, after which pull requests are made for merging with `devel`.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* A merge from `devel` to the `master` branch invokes a new release.

## Contributors
* [Ìván Domenzain](https://www.chalmers.se/en/staff/Pages/ivand.aspx) ([@IVANDOMENZAIN](https://github.com/IVANDOMENZAIN)), Chalmers University of Technology, Göteborg, Sweden
* [Linnea Österberg](https://www.chalmers.se/en/search/Pages/default.aspx?q=linnea+%c3%b6sterberg) ([@linoste](https://github.com/linoste)), Chalmers University of Technology, Göteborg, Sweden
