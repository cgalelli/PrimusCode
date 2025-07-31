# PrimusCode
**A Simulation Pipeline for Cosmic Ray Paleo-Detectors** 

This repository contains the Python-based analysis framework for the **PRImuS** (Paleo-astroparticles Reconstructed with the Interactions of MUons in Stone) project. The code is designed to perform phenomenological studies on paleo-detectors, specifically by simulating the expected rate of nuclear recoil tracks induced by cosmic ray muons in natural minerals. The **primary scientific application** of this code is detailed in our paper, available on **arXiv**: [https://arxiv.org/abs/2405.04908](https://arxiv.org/abs/2405.04908)

## The Analysis Workflow
The analysis pipeline is structured around a single Jupyter Notebook (`Primus_analysis.ipynb`) that interacts with a core utility module (`mineral_utils.py`). All parameters for a given analysis run are defined in a central configuration dictionary at the beginning of the notebook. This includes:
 - The mineral to be analyzed (e.g., Halite, Olivine).
 - The geological history of the sample (age, exposure time, deposition rate).
 - The astrophysical scenarios for the cosmic ray flux (e.g., standard flux, enhanced flux from a nearby supernova).

To run this analysis pipeline, you will need a Python environment with the following packages installed. You can install them using `pip`:
```bash
pip install numpy scipy matplotlib mendeleev
```
To run the analysis, clone the Repository:
```bash
git clone https://github.com/cgalelli/PrimusCode.git
cd PrimusCode
```
## Citation

If you use this work or the associated data in your research, please cite our paper:
```yaml
@article{Caccianiga_2024,
    author        = {Caccianiga, Lorenzo and Apollonio, Lorenzo and Mariani, Federico Maria and Magnani, Paolo and Galelli, Claudio and Veutro, Alessandro},
    title         = {{Sedimentary rocks from Mediterranean drought in the Messinian age as a probe of the past cosmic ray flux}},
    journal       = {Phys. Rev. D},
    volume        = {110},
    issue         = {12},
    pages         = {L121301},
    year          = {2024},
    doi           = {10.1103/PhysRevD.110.L121301}
}
```

## Contact

For questions, please contact: Claudio Galelli – claudio.galelli@mi.infn.it
