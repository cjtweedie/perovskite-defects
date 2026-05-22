# perovskite-defects
Metal-halide perovskites with chemical formulas ABX3 are important energy materials, but to improve their efficiancy and stability we need to understand the microscopic behaviour. This can be done through atomistic modelling from first-principles (e.g. DFT) or through machine-learned forced fields (e.g. MLFF-MD). For now, everything will focus on interfacing with DFT calculations in VASP, so input/output files are expected to be compatible. Dependency on ASE, SNB, ICET python packages.

The primary function of these scripts is for vacancy and interstitial generation, especially for intrinsic halide defects. Finding the most favourable interstitial configuration for atomistic modelling is difficult, so a systematic potential energy surface scan with ionic relaxation is used here to consider all possible metastable and stable states and compare their relative energies. Electronic/vibrational structure can also be analysed from DFT calculations of these defect systems. 

There will be scripts to set up transition-state calculations for defect-mediated ion migration, either with DFT+NEB or MLFF+NEB. I might try to find a way to interface with Rohan's code so it's one big workflow.

There are also scripts for modelling mixed-composition perovskites within the best finite-size supercell approximation, using a stochastic approach to the special quasirandom structures (SQS) method (LINK TO PAPER?) enabled by the ICET package (LINK TO THIS TOO). Mixing of A or B (metal/organic) cations or X (halide) ions is available, and simple binary or multi-component alloys can be created in one or multiple sublattices of the structure. Input parameters for the cluster expansion  can be left as default or tuned to improve SQS accuracy/reproducibility. 

Probably some more stuff here too when I get to it.
