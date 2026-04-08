# PrimusCode
**A Simulation Pipeline for Cosmic Ray Paleo-Detectors** 

This repository contains the Python-based analysis framework for the **PRImuS** (Paleo-astroparticles Reconstructed with the Interactions of MUons in Stone) project. The code is designed to perform phenomenological studies on paleo-detectors, specifically by simulating the expected rate of nuclear recoil tracks induced by cosmic ray muons in natural minerals. A first **scientific application** of this code is detailed in our paper, available on **arXiv**: [https://arxiv.org/abs/2405.04908](https://arxiv.org/abs/2405.04908)

## The Analysis Workflow
The analysis pipeline is structured around Jupyter Notebooks (such as `CDP_analysis_new.ipynb`) that interacts with a core utility module (`mineral_utils.py`). All parameters for a given analysis run are defined in a central configuration dictionary, mostly imported using `yaml` at the beginning of the notebook. This includes:
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
@article{Galelli:2025gss,
    author = "Galelli, Claudio and Caccianiga, Lorenzo and Apollonio, Lorenzo and Magnani, Paolo and Breton, Vincent",
    title = "{A volcanic chronosequence as a time-resolved paleo-detector array to study the cosmic-ray flux in the Late Pleistocene and Holocene}",
    eprint = "2510.23126",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.HE",
    month = "10",
    year = "2025"
}
```

Or the associated proceedings from ICRC 2025
```yaml
@article{Galelli:2025,
    author        = {Galelli, Claudio and Caccianiga, Lorenzo and Apollonio, Lorenzo and Magnani, Paolo},
    title         = {{Probing Ancient Cosmic Ray Flux with Paleo-Detectors and the Launch of the PRImuS Project}},
    journal       = {Proceedings of the 39th International Cosmic Ray Conference, ICRC 2025},
    year          = {2025},
}
```

## Contact

For questions, please contact: Claudio Galelli – claudio.galelli@mi.infn.it
