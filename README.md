# Search for ALP in associated vector boson production at LHC

## Production in Madgraph 
Process $pp \to ALP/H (\gamma \gamma) W (\ell \nu)$ is generated in Madgraph for 14 TeV CM energy and later for hadronization Pythia8 and collider efficiencies are imposed using Delphes (CMS card). Similarly, processes including $Z$ boson are also produced where $Z \to \ell^+ \ell^-$. 

LHE files for ALP and Higgs for pp > ALP/Higgs Z process at 14 TeV CM energy are here for 125 & 25 GeV mass point  https://www.dropbox.com/scl/fo/cdlnwpd9sl4dxnhvcw2fh/ANgcBYqkYiWP2-rzO9eesc0?rlkey=15syagbqxxz6kft4mibadwxcs&st=ppydn3el&dl=0

## Delphes analyzer
To run the Delphes analyzer for the Z-related process with the root files of the individual processes, keep the ```ALP_AP_Z.C``` in the ```Delphes.X.X``` folder and run as ```root -l ALP_AP_Z.C'("/path_to_file/XX.root")'```  
## XGBoost Analysis
XGBoost analysis for ALP vs Higgs at the 125 GeV mass point is in ```XGBoost_ALP_AP.ipynb```. 
