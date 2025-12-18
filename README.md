This repository contains computational supporting information for the manuscript "Importance of Nonadiabatic Effects in Modeling Absorption Spectra of Indole and Cyanoindole Derivatives in Solution" by D. Bashirova, J. Stein, K. Siegfried, and T. J. Zuehlsdorff
<p align="center">
  <img src="https://github.com/user-attachments/assets/afe1f4e7-bb78-401d-87c1-4a09ccd77682" width='300' />
</p>

```bash
.
├── 2CNI
│   ├── GCT  # contains raw absorption spectra and spectral densities
│   ├── GNCT # contains raw absorption spectra and spectral densities
│   │   ├── s1
│   │   └── s2
│   ├── MD # contains topology and parameter files, as well as excitation energies along the QM/MM simulation 
│   └── T-TEDOPA
│       ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│       └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
├── 4CNI
│   ├── GCT  # contains raw absorption spectra and spectral densities
│   ├── GNCT # contains raw absorption spectra and spectral densities
│   │   ├── s1
│   │   └── s2
│   ├── MD # contains topology and parameter files, as well as excitation energies along the QM/MM simulation
│   └── T-TEDOPA
│       ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│       └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
├── 6CNI
│   ├── GCT  # contains raw absorption spectra and spectral densities
│   ├── GNCT # contains raw absorption spectra and spectral densities
│   │   ├── s1
│   │   └── s2
│   ├── MD # contains topology and parameter files, as well as excitation energies along the QM/MM simulation
│   └── T-TEDOPA
│       ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│       └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
├── 7CNI
│   ├── GCT  # contains raw absorption spectra and spectral densities
│   ├── GNCT # contains raw absorption spectra and spectral densities
│   │   ├── s1
│   │   └── s2
│   ├── MD # contains topology and parameter files, as well as excitation energies along the QM/MM simulation
│   └── T-TEDOPA
│       ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│       └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
├── Indole
│   ├── Vacuum
│   │   ├── AIMD # contains an input geometry and excitation energies along the AIMD simulation
│   │   ├── GCT  # contains raw absorption spectra and spectral densities
│   │   ├── GNCT # contains raw absorption spectra and spectral densities
│   │   │   ├── s1
│   │   │   └── s2
│   │   └── T-TEDOPA
│   │       ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│   │       └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
│   └── Water
│       ├── GCT  # contains raw absorption spectra and spectral densities
│       ├── GNCT # contains raw absorption spectra and spectral densities
│       │   ├── s1
│       │   └── s2
│       ├── MD # contains topology and parameter files, as well as excitation energies along the QM/MM simulation
│       └── T-TEDOPA
│           ├── Diabatic_data # contains diabatic energies, dipole moments, and couplings, as well as diabatic spectral densities
│           └── Quantum_Dynamics # contains the input file, chain coefficients, resulting absorption spectrum, and population dynamics
├── Input_files # contains input files for AIMD, MM, and QM/MM simulations; TDDFT excitation energies; and GCT/GNCT absorption spectra and adiabatic spectral density calculations
└── Python_scripts # contains diabatization, diabatic spectral density generation, and chain mapping scripts

