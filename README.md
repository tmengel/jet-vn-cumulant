# Jet $v_{n}$ Cumulant Method

This framework creates toy-model simulations of 200 GeV Au+Au heavy ion collisions and calcualtes jet $v_{n}$ through mixed cumulants.

## Description

The toy-model created in simulation takes a stand alone version of TennGen to create 200 GeV Au+Au backgrounds with a set event plane and particle $v_{N}$. The background
can be made for 0-10%, 10-20%, 20-40% and 40-60% central collisions with a maximum psuedo-rapidity of $|\eta| < 1.1$. The background only consists of $\pi^{\pm}$, $k^{\pm}$, $p$, and $\bar{p}$. 
The signal jets are simualted with PYTHIA 8 200 GeV $p$+$p$ collisions. All final state particles are saved from the PYTHIA event. The FastJet3 is used to find anti-$k_{T}$ jets
for $R = 0.4$ jets. The leading jet is then aligned to match a $\frac{dN_{jet}}{d\phi}$ distribution where the event-plane is defined by the TennGen event plane and the jet $v_{n}$ is set by the user.
The entire Pythia event is then rotated and merged with the TennGen background. 

The analyis macro takes in the mixed, re-aligned event and finds anti-$k_{T}$ jets with $R=0.4$ and a $p_{T}^{min} = 10$ GeV.

## Getting Started

### Dependencies

* PYTHIA 8.3.X
* ROOT 6.*X
* FastJet 3.X

### Build

* Update the Makefile.Inc to you own local installations of the required dependencies.
* cd to any Macro directory and type <make> 

### Executing simulation macros
 
* ./GenerateTennGenAuAu <nevents> <centrality> <eta range> <pt bin to match with pythia>
  - Creates root files for background in /simulation/root-files/AuAu/<centrality>
* ./GeneratePythiaPP <Collision Energy> <nevents> <ptbin>
  - Creates root files for signal in /simulation/root-files/PP
* ./RealignJets <input pythia file> <input tenngen file> <prefix> <jet resolution parameter> <nEvents>
  - Roatates pythia event and merges with TennGen background, outputs to <prefix>/200GeV_MixedEvents_ptbin<ptbin of input file>_RealignJets_R<jet param>.root
* ./MergeAllPtBins <input_dir> <n_pt_bins> <jetparam> <output_file>
  - Merges all pt bins into one file in /simulation/root-files/MixedEvents/<output_file>
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

ex. Tanner Mengel
ex. [tmengel@vols.utk.edu]

## Version History

* 0.1
    * Creation
