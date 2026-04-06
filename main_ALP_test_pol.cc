////https://feynrules.irmp.ucl.ac.be/wiki/ALPsEFT


#include "Pythia8/Pythia.h"
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif
// ROOT, for histogramming.
//#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"

using namespace Pythia8;

// ROOT, for saving file.
#include "TFile.h"
#include "TLorentzVector.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include <fastjet/tools/MassDropTagger.hh>
//#include "fastjet/contrib/ModifiedMassDropTagger.hh"
//#include "fastjet/contrib/SoftDrop.hh"
using namespace fastjet;

////////Function to calculate reco nu 
TLorentzVector reconstructNeutrino(const TLorentzVector& lep, double met_px, double met_py) {

    const double MW = 80.4;

    double px_l = lep.Px();
    double py_l = lep.Py();
    double pz_l = lep.Pz();
    double E_l  = lep.E();

    double pt_l2 = px_l*px_l + py_l*py_l;

    double Lambda = MW*MW/2.0 + px_l*met_px + py_l*met_py;

    double discriminant = Lambda*Lambda - pt_l2*(met_px*met_px + met_py*met_py);

    double pz_nu;

    if (discriminant >= 0) {
        double sqrtD = sqrt(discriminant);

        double pz1 = (Lambda * pz_l + E_l * sqrtD) / pt_l2;
        double pz2 = (Lambda * pz_l - E_l * sqrtD) / pt_l2;

        // choose smaller |pz| or you can choose the one which W mass closest to MW value 
        pz_nu = (fabs(pz1) < fabs(pz2)) ? pz1 : pz2;
    } else {
        //real part only
        pz_nu = (Lambda * pz_l) / pt_l2;
    }

    TLorentzVector nu;
    double Enu = sqrt(met_px*met_px + met_py*met_py + pz_nu*pz_nu);
    nu.SetPxPyPzE(met_px, met_py, pz_nu, Enu);

    return nu;
}


int main() {

  //Pythia8::Pythia8ToHepMC topHepMC("/Users/amandip/delphes-master/hepmcout_gamgam_A10.dat");
  //Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Next:numberCount = 10000");
  
  // Read in commands from external file.
  //pythia.readFile("main1000.cmnd");
  //pythia.readString("PDF:pSet = 14"); 
  //pythia.readString("Tune:pp=21");
  // Pick new random number seed for each run, based on clock.
  //pythia.readString("Random:setSeed = on");
  //pythia.readString("Random:seed = 0");
  //pythia.readString("PhaseSpace:mHatMin = 500.");
  //pythia.readString("PhaseSpace::pTHatMin = 500.");
  // Extract settings to be used in the main program.
  //int nEvent   = pythia.mode("Main:numberOfEvents");
///home/amandeep/MG5_aMC_v3_5_7/bin/TEST_VBF_WW/Events/run_01

  // Read from LHE file 
  pythia.readString("Beams:frameType = 4");
  
  //pythia.readString("Beams:LHEF =/data/SOM_Test/ALP_test/AP_Higgs_wm_lep_gam/Events/run_02/unweighted_events.lhe");
  //pythia.readString("Beams:LHEF =/data/SOM_Test/ALP_test/AP_Higgs_wp_lep_gam/Events/run_01/unweighted_events.lhe");
  
  pythia.readString("Beams:LHEF =/data/SOM_Test/ALP_test/AP_ALP_Wp_lep_gam/Events/run_02/unweighted_events.lhe"); // W+
  //pythia.readString("Beams:LHEF =/data/SOM_Test/ALP_test/AP_ALP_Wm_lep_gam/Events/run_02/unweighted_events.lhe"); // W-

  pythia.init();
  int nEvent = 100000;
  //int nEvent   = pythia.mode("Main:numberOfEvents");

  TH1F* missE_all = new TH1F("missing Et", "Missing Et", 100,0.0,1000);
  TH1F* hist_neutrino_pt = new TH1F("neutrino_pt", "neutrino_pt", 100,0.0,1000);


  //TH1F *hist_= new TH1F("delta R of two gamma from ALP","delta R of two gamma from ALP", 40, 0.0, 1.);
  TH1F *hist_W_pt = new TH1F("pt_W","pt_W", 100, 0.0, 1000.);
  TH1F *hist_W_pt1 = new TH1F("pt_W1","pt_W1", 100, 0.0, 1000.);
  TH1F *hist_W_pt2 = new TH1F("pt_W2","pt_W2", 100, 0.0, 1000.);
  TH1F *hist_W_pt3 = new TH1F("pt_W3","pt_W3", 100, 0.0, 1000.);
  TH1F *hist_dphi_Wrecopt_digam = new TH1F("dphi_Wrecopt_digam","dphi_Wrecopt_digam", 100, -6.0, 6.);
  TH1F *hist_deta_Wrecopt_digam = new TH1F("deta_Wrecopt_digam","deta_Wrecopt_digam", 100, -6, 6.);
  TH1F *hist_LP = new TH1F("LP","LP", 100, -2, 2.);
  TH1F *hist_ratio_ptl_ptW = new TH1F("ratio_ptl_ptW","ratio_ptl_ptW", 100, 0, 1);
  TH1F *hist_ratio_ptl_ptdigam = new TH1F("ratio_ptl_ptdigam","ratio_ptl_ptdigam", 100, 0, 1);





  TH1F *hist_dphi_Wrecopt_W = new TH1F("dphi_Wrecopt_W","dphi_Wrecopt_W", 100, -6.0, 6.);
  TH2F* h2_ptlep_vs_Wreco = new TH2F("ptlep_vs_Wreco", "p_{T}^{l} vs reco p_{T}^{W}; reco p_{T}^{W} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  TH2F* h2_ptlep_vs_digam = new TH2F("ptlep_vs_digam", "p_{T}^{l} vs p_{T}^{2#gamma}; p_{T}^{2#gamma} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  

  TH1F *hist_W_eta = new TH1F("eta_W","eta_W", 100, -6.0, 6.);
  TH1F *hist_W_mass = new TH1F("m_W","m_W", 100, 0, 100);
  TH1F *hist_W_phi = new TH1F("phi_W","phi_W", 100, -6.0, 6.);
  TH1F* h_cos_theta_star = new TH1F("cos_theta_star", "cos_theta_star", 20, -1.0, 1.0);  
  TH1F* h_cos_theta_star_cut1 = new TH1F("cos_theta_star_cut1", "cos_theta_star1", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_cut2 = new TH1F("cos_theta_star_cut2", "cos_theta_star2", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_cut3 = new TH1F("cos_theta_star_cut3", "cos_theta_star3", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_cut4 = new TH1F("cos_theta_star_cut4", "cos_theta_star4", 20, -1.0, 1.0); 
  
  TH1F* h_phi_star_cut = new TH1F("phi_star_cut", "phi_star_cut", 20, 0., 3.5); 
  TH1F* h_phi_star_cut1 = new TH1F("phi_star_cut1", "phi_star_cut1", 20, 0., 3.5); 
  TH1F* h_phi_star_cut2 = new TH1F("phi_star_cut2", "phi_star_cut2", 20, 0., 3.5);
  TH1F* h_phi_star_cut3 = new TH1F("phi_star_cut3", "phi_star_cut3", 20, 0., 3.5);
  TH1F* h_phi_star_cut4 = new TH1F("phi_star_cut4", "phi_star_cut4", 20, 0., 3.5);


  TH1F* h_cos_theta1 = new TH1F("cos_theta1", "cos_theta1", 20, -1.0, 1.0);  


  TH1F* h_theta_star = new TH1F("theta_star","theta_star", 20, 0.0, 180.0);
  TH1F *hist_cut_flow = new TH1F("cut flow", "cut flow", 6, 0, 6);

  TH2F* h2_cos_vs_ptW = new TH2F("cos_vs_ptW", "cos#theta* vs p_{T}^{W};p_{T}^{W} [GeV];cos#theta*", 50, 0, 250, 50, -1, 1);
  TH2F* h2_cos_vs_ptlep = new TH2F("cos_vs_ptlep", "cos#theta* vs p_{T}^{lep};p_{T}^{lep} [GeV];cos#theta*", 50, 0, 100, 50, -1, 1);

  TH1F *hist_dphi_lep_gamma1 = new TH1F("dphi_lep_gamma1","dphi_lep_gamma1", 20, 0.0 , 3.5);
  TH1F *hist_dphi_lep_gamma2 = new TH1F("dphi_lep_gamma2","dphi_lep_gamma2", 20, 0.0 , 3.5);
  TH1F *hist_dR_lep_gamma1 = new TH1F("dR_lep_gamma1","dR_lep_gamma1", 20, 0.0 , 3.5);



  TH1F *hist_lep_pt = new TH1F("lep_pt","lep_pt", 25, 0.0, 250.);
  TH1F *hist_dphi_lep_MET = new TH1F("dphi_lep_MET","dphi_lep_MET", 100, -6 , 6);
  TH1F *hist_sig_lep_rap = new TH1F("sig_lep_rap","sig_lep_rap", 20, -4.5 , 4.5);



  //////////////////Varibles ///////
   //TFile *outfile = new TFile("/data/SOM_Test/ALP_test/OutputFile_Higgs_Wp_125.root","RECREATE");

  TFile *outfile = new TFile("/data/SOM_Test/ALP_test/OutputFile_ALP_Wp_125.root","RECREATE");
  TTree *tree=new TTree("Tree","Signal");
  Float_t  W_reco_pt = 0., Lep_pt=0., MET_et = 0.0,  Cos_theta_star =0.0, Dphi_Lep_MET=0., Dphi_Wrecp_Diphoton=0.,Deta_Wrecp_Diphoton=0., Phi_star=0.;
  Float_t  LP=0.;

  tree->Branch("W_reco_pt",&W_reco_pt,"W_reco_pt/F");
  tree->Branch("Lep_pt",&Lep_pt,"Lep_pt/F");
  tree->Branch("MET_et",&MET_et,"MET_et/F");
  tree->Branch("Cos_theta_star",&Cos_theta_star,"Cos_theta_star/F");
  tree->Branch("Dphi_Lep_MET",&Dphi_Lep_MET,"Dphi_Lep_MET/F");
  tree->Branch("Cos_theta_star",&Cos_theta_star,"Cos_theta_star/F");
  tree->Branch("Dphi_Wrecp_Diphoton",&Dphi_Wrecp_Diphoton,"Dphi_Wrecp_Diphoton/F");
  tree->Branch("Deta_Wrecp_Diphoton",&Deta_Wrecp_Diphoton,"Deta_Wrecp_Diphoton/F");
  tree->Branch("Phi_star",&Phi_star,"Phi_star/F");
  tree->Branch("Dphi_Wrecp_Diphoton",&Dphi_Wrecp_Diphoton,"Dphi_Wrecp_Diphoton/F");
  tree->Branch("LP",&LP,"LP/F");



  ////////////////////////////



int n1 = 0 , n2 = 0 , n3 = 0, n4 =0 , n5 = 0 , n6 =0; 

const int    nPtBins = 20;
const double ptMin  = 0.0;
const double ptMax  = 500.0;
const double ptBinW = (ptMax - ptMin) / nPtBins;


const int    netaBins = 16;
const double etaMin  = -4;
const double etaMax  = +4;
const double etaBinW = (etaMax - etaMin) / netaBins;

std::vector<TH1F*> h_costheta_wpt;
for (int i = 0; i < nPtBins; ++i) {

  double lo = ptMin + i * ptBinW;
  double hi = lo + ptBinW;
  h_costheta_wpt.push_back(
    new TH1F(
      Form("cos_theta_star_wptbin_%d", i),
      Form("cos#theta^{*}, %.0f < p_{T}^{W} < %.0f;cos#theta^{*};Events", lo, hi),
      20, -1.0, 1.0
    )
  );
}
std::vector<TH1F*> h_costheta_leta;
for (int i = 0; i < netaBins; ++i) {

  double lo = etaMin + i * etaBinW;
  double hi = lo + etaBinW;
 h_costheta_leta.push_back(
    new TH1F(
      Form("cos_theta_star_etabin_%d", i),
      Form("cos#theta^{*}, %.0f < #eta^{l} < %.0f;cos#theta^{*};Events", lo, hi),
      20, -1.0, 1.0
    )
  );
}

std::vector<TH1F*> h_costheta_lpt;
for (int i = 0; i < nPtBins; ++i) {

  double lo = ptMin + i * ptBinW;
  double hi = lo + ptBinW;
  h_costheta_lpt.push_back(
    new TH1F(
      Form("cos_theta_star_lptbin_%d", i),
      Form("cos#theta^{*}, %.0f < p_{T}^{l} < %.0f;cos#theta^{*};Events", lo, hi),
      20, -1.0, 1.0
    )
  );
}
/////////////////////Asymmetry with W pt /////////////
const int nAFBbins = 4;
double wptBins[nAFBbins+1] = {30, 50, 80, 120, 700};
double Nf[nAFBbins] = {0.0};  // forward
double Nb[nAFBbins] = {0.0};  // backward

//////////////////
  //--------------------------------------------start pythia loop-------------------------------------------------
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

   //if (iEvent == 60) pythia.event.list() ;
    //if (iEvent >2) continue; 

// vector<fastjet::PseudoJet> ALP_gamma;
// vector<fastjet::PseudoJet> ALP_others;

// vector<fastjet::PseudoJet> fjInputs, PARTICLES_nu, VBF_JETS, fjInputs_1;
// Keep track of missing ET
  Vec4 missingETvec;
  double px_met = 0, py_met = 0;
//topHepMC.writeNextEvent( pythia );



// Find W boson and its decay products
        int W_idx = -1;
        int lepton_idx = -1;
        int neutrino_idx = -1;
        int ALP_idx = -1;
        int light_dau1 = -1;
        int light_dau2 = -1;
        int ALP_gam1 = -1;
        int ALP_gam2 = -1;

//--------------------------------------------start pythia event looop-------------------------------------------------
for (int i = 0; i < event.size(); ++i){
    auto &p = event[i];
    fastjet::PseudoJet pj(p.px(),p.py(),p.pz(),p.e());

//if (p.idAbs() == 24) cout<<p.idAbs()<<"  "<<p.status()<<"   "<<iEvent<<endl; 
//if (p.idAbs() == 9000005 || p.idAbs() == 25) cout<<p.idAbs()<<"  "<<p.status()<<"   "<<iEvent<<endl;

if (p.idAbs() == 24 && abs(p.status()) == 62) {  // W boson (after FSR)
                W_idx = i;
}


if ( (p.idAbs() == 25 || p.idAbs() == 9000005 || p.idAbs() == 36) && abs(p.status()) == 62) {  // ALP
                ALP_idx= i;
}


if (p.idAbs() == 24 && abs(p.status()) == 62) {  // W boson (after FSR)
                W_idx = i;
                // Check daughters
                int d1 = p.daughter1();
                int d2 = p.daughter2();
                
                if (d1 > 0 && d2 > 0) {
                    Particle& daughter1 = event[d1];
                    Particle& daughter2 = event[d2];
                    
                    // Check for leptonic decay (e, mu, tau)
                    if (daughter1.idAbs() == 11 || daughter1.idAbs() == 13 || 
                        daughter1.idAbs() == 15) {
                        lepton_idx = d1;
                        neutrino_idx = d2;

                          //cout<<d1<<"  "<<d2<<" "<<iEvent<<endl;
                   
                    } else if (daughter2.idAbs() == 11 || daughter2.idAbs() == 13 || 
                               daughter2.idAbs() == 15) {
                        lepton_idx = d2;
                        neutrino_idx = d1;
                    
                    }
                  // If the decay of w+ to j j 
                   else if ((daughter1.idAbs() >= 1 && daughter1.idAbs() <= 4) || (daughter2.idAbs() >= 1 && daughter2.idAbs() <= 4))
                   {
                    light_dau1 = d1;
                    light_dau2 = d2;
                   }


                }
              }



if ( (p.idAbs() == 25 || p.idAbs() == 9000005) ) {
  
   auto d1 = pythia.event[p.daughter1()];
   auto d2 = pythia.event[p.daughter2()]; 
 
    //cout<<abs(d1.eta() - d2.eta()) <<"  "<<d1.idAbs()<<"  "<<d2.idAbs()<<"  "<<iEvent<<endl;
    //if (d1.idAbs() !=22 && d1.idAbs() !=22 )cout<<d1.id()<<"  "<<d2.id()<<"  "<<iEvent<<endl;

  bool daughters_valid = true;
    
    if (p.idAbs() == 25 && (d1.idAbs() == 25 || d2.idAbs() == 25)) {
        daughters_valid = false;
    }
    if (p.idAbs() == 9000005 && (d1.idAbs() == 9000005 || d2.idAbs() == 9000005)) {
        daughters_valid = false;
    }

    if(daughters_valid == false) continue;

   if (d1.idAbs() == 22) {
            //cout<<d1.id()<<"   "<<d1.status()<<"   " <<p.id()<<"  "<<iEvent<<endl;
            ALP_gam1 = p.daughter1();
        }

   if (d2.idAbs() == 22){ 
            //cout<<d2.id()<<"   "<<d2.status()<<"   " <<p.id()<<"  "<<iEvent<<endl;
            ALP_gam2 = p.daughter2();
      }
  
   //cout<<d1.id()<<"   "<<d2.id()<<"   " <<p.id()<<"  "<<iEvent<<endl;
   // ALP_others.push_back(fastjet::PseudoJet (d1.px(),d1.py(),d1.pz(),d1.e() ));
   // ALP_others.push_back(fastjet::PseudoJet (d2.px(),d2.py(),d2.pz(),d2.e() ));
}

//Final state particles
if (! p.isFinal() ) continue ; // Final State only 
//cout<<p.status()<<endl;
///////////Detector like MET ; MET = - visible sum 
// No neutrinos
if (p.idAbs() == 12 || p.idAbs() == 14 || p.idAbs() == 16)  continue;

px_met -= pythia.event[i].px();
py_met -= pythia.event[i].py();

// construct the particle for missing ET
  //PARTICLES_nu.push_back(fastjet::PseudoJet(p.px(),p.py(),p.pz(),p.e()));

} // End pythia event loop
  
n1++;
//// Missing Energy
double met_et=0.0;
met_et = sqrt(px_met*px_met + py_met*py_met);
//missE_all->Fill(met_et); ///var

//cout<<px_met<<"   "<<py_met<<"   "<<sqrt(px_met*px_met + py_met*py_met)<<endl;

//cout<<ALP_idx<<"    "<<ALP_gam1<<"  "<<ALP_gam2<<"   "<<W_idx<<" "<<lepton_idx<<"  "<<neutrino_idx<<"  "<<light_dau1<<"   "<<light_dau2<<"   "<<iEvent<<endl;


if(W_idx >= 0 &&  lepton_idx >= 0 && neutrino_idx >= 0 && ALP_gam1 >=0 && ALP_gam2 >= 0){

n2++;

Particle& W = event[W_idx];
Particle& ALP = event[ALP_idx];
Particle& lepton = event[lepton_idx];
Particle& neutrino = event[neutrino_idx];
Particle& gamma1 = event[ALP_gam1];
Particle& gamma2 = event[ALP_gam2];

TLorentzVector MET_vec(px_met, py_met, 0.0, met_et);
TLorentzVector W_vec(W.px(), W.py(), W.pz(), W.e());
TLorentzVector lep_vec(lepton.px(), lepton.py(), lepton.pz(), lepton.e());
TLorentzVector lep_vec_org(lepton.px(), lepton.py(), lepton.pz(), lepton.e());

TLorentzVector neu_vec(neutrino.px(), neutrino.py(), neutrino.pz(), neutrino.e());
TLorentzVector gam1_vec(gamma1.px(), gamma1.py(), gamma1.pz(), gamma1.e());
TLorentzVector gam2_vec(gamma2.px(), gamma2.py(), gamma2.pz(), gamma2.e());
TLorentzVector diphoton = gam1_vec + gam2_vec;
// reconstruct neutrino from met 
TLorentzVector nu_reco = reconstructNeutrino(lep_vec_org, px_met, py_met);
// reconstructed W
TLorentzVector W_reco = lep_vec_org + nu_reco;

//////////W reco
double cos_theta_star =-5, phi_star = -5;

if (lepton.pT() > 20 && abs(lepton.eta()) < 2.5  && gamma1.pT() > 10 && gamma2.pT() > 10 && abs(gamma1.eta()) < 2.5 && abs(gamma2.eta()) < 2.5  && met_et > 20){ 
n3++;  
missE_all->Fill(met_et); ///var
W_reco_pt = 0., Lep_pt = 0., MET_et=0., Dphi_Wrecp_Diphoton=0., Deta_Wrecp_Diphoton=0.;
W_reco_pt = W_reco.Pt();
Lep_pt = lepton.pT();
MET_et = met_et;
Dphi_Wrecp_Diphoton=abs(W_reco.DeltaPhi(diphoton));
Deta_Wrecp_Diphoton=abs(W_reco.Eta() - diphoton.Eta());
//Gen W 
hist_W_pt->Fill(W.pT());
hist_W_eta->Fill(W.eta());
hist_W_mass->Fill(W.m());
hist_W_phi->Fill(W.phi());

hist_W_pt1->Fill(diphoton.Pt());
hist_W_pt2->Fill((MET_vec + lep_vec_org).Pt()); // Vec sum of MET & Lep 
hist_W_pt3->Fill(W_reco.Pt()); /////// var 
hist_dphi_Wrecopt_digam->Fill(Dphi_Wrecp_Diphoton);
hist_deta_Wrecopt_digam->Fill(Deta_Wrecp_Diphoton);
hist_dphi_Wrecopt_W->Fill(W_reco.DeltaPhi(W_vec));
h2_ptlep_vs_Wreco->Fill(lepton.pT(), W_reco.Pt());
h2_ptlep_vs_digam->Fill(lepton.pT(), diphoton.Pt());
//cout<<W.pT()<<"   "<<(MET_vec + lep_vec).Pt()<<"    "<<W_reco.Pt()<<"   "<<diphoton.Pt()<<endl;
//cout<<W_reco.Px()<<"   "<<nu_reco.Px()<<"   "<<lep_vec.Px()<<endl;


/////////////////////
/* Using gen W for now ; replace reco W later*/
TVector3 boostW = W_vec.BoostVector();
lep_vec.Boost(-boostW);

TVector3 lep_dir_Wrest = lep_vec.Vect().Unit();   // lepton in W rest frame
TVector3 W_dir_lab     = W_vec.Vect().Unit();     // W direction in lab

cos_theta_star= lep_dir_Wrest.Dot(W_dir_lab);
//h_cos_theta_star->Fill(cos_theta_star);

////////////Costeta 1 
    // CM frame of final state objetcs
    TLorentzVector ALL_final = W_reco + diphoton;
    TVector3 boostCM = ALL_final.BoostVector();
    // To get the W dir in the CM frame
    TLorentzVector W_reco_CM = W_reco ;
    W_reco_CM.Boost(-boostCM);
    // Dir of W in the CM frame
    TVector3 W_dir_CM = W_reco_CM.Vect().Unit();
    double cos_theta1 = W_dir_CM.CosTheta();  // wrt z-axis
    // Dir of Higgs in the CM frame
    TLorentzVector diphoton_CM = diphoton;
    diphoton_CM.Boost(-boostCM);
    TVector3 H_dir_CM = diphoton_CM.Vect().Unit();

    // check 
    //cout<<" In CM this should be -1 => "<<W_dir_CM.Dot(H_dir_CM)<<endl;
    h_cos_theta1->Fill(cos_theta1);

///////////////////// Polarization Fraction ////////////
double y_lep = -5, pt_W=-5, pt_lep = -5; 
y_lep = lepton.y(); // rapidity 
pt_W = W.pT();
pt_lep = lepton.pT();

//if (lepton.pT() > 20 && abs(lepton.eta()) < 2.5  && gamma1.pT() > 10 && gamma2.pT() > 10 && abs(gamma1.eta()) < 2.5 && abs(gamma2.eta()) < 2.5 && met_et > 20){ 
//if(diphoton.Pt() > 100){
Cos_theta_star=0.;
Cos_theta_star= cos_theta_star;

h_cos_theta_star->Fill(cos_theta_star);

int ipt_W = int((pt_W - ptMin) / ptBinW);
if (ipt_W >= 0 && ipt_W < nPtBins) {
  h_costheta_wpt[ipt_W]->Fill(cos_theta_star);
}
int ieta = int((y_lep - etaMin) / etaBinW);
if (ieta >= 0 && ieta < netaBins) {
    h_costheta_leta[ieta]->Fill(cos_theta_star);
}

int ipt_lep = int((pt_lep - ptMin) / ptBinW);
if (ipt_lep >= 0 && ipt_lep < nPtBins) {
  h_costheta_lpt[ipt_lep]->Fill(cos_theta_star);
}
////////////////////////////AFB///////////////

int ibin = -1;
double pt_W_reco = W_reco.Pt();
for (int i = 0; i < nAFBbins; ++i) {
    if (pt_W_reco >= wptBins[i] && pt_W_reco < wptBins[i+1]) {
        ibin = i;
        break;
    }
}

if (ibin >= 0) {
    if (cos_theta_star > 0) {
        Nf[ibin] += 1.0;
    } else {
        Nb[ibin] += 1.0;
    }
}

//////////////////////////////

//phi_star azimuthal angle of the charged lepton in the W rest frame
//measured with respect to the production plane
TVector3 beam_dir(0.0, 0.0, 1.0);   // proton beam direction
TVector3 y_hat = beam_dir.Cross(W_dir_lab).Unit();
TVector3 x_hat = y_hat.Cross(W_dir_lab).Unit();

phi_star = abs(atan2(lep_dir_Wrest.Dot(y_hat),lep_dir_Wrest.Dot(x_hat)));
h_phi_star_cut->Fill(phi_star);
//////////
//if (lepton.pT() >= 20 && abs(lepton.eta()) <= 4.5 && neutrino.pT() >= 20) {

h2_cos_vs_ptW->Fill(W.pT(), cos_theta_star);
h2_cos_vs_ptlep->Fill(lepton.pT(), cos_theta_star);

hist_lep_pt->Fill(lepton.pT());
hist_neutrino_pt->Fill(neutrino.pT());
hist_sig_lep_rap->Fill(lepton.eta());

hist_dphi_lep_MET->Fill(abs(lep_vec_org.DeltaPhi(nu_reco)));

hist_dphi_lep_gamma1->Fill(lep_vec_org.DeltaPhi(gam1_vec));
hist_dphi_lep_gamma2->Fill(lep_vec_org.DeltaPhi(gam2_vec));
hist_dR_lep_gamma1->Fill(lep_vec_org.DeltaR(gam2_vec));
//}
Dphi_Lep_MET=0., Phi_star=0.;
Dphi_Lep_MET=abs(lep_vec_org.DeltaPhi(nu_reco));
Phi_star = phi_star;
////////////L_p
LP=0.;
TVector3 pT_3lep = lep_vec_org.Vect();
TVector3 pT_3W   = W_reco.Vect();
LP= pT_3lep.Dot(pT_3W)/pT_3W.Mag2();
hist_LP->Fill(LP);


hist_ratio_ptl_ptW->Fill(lepton.pT()/W_reco.Pt());
hist_ratio_ptl_ptdigam->Fill(lepton.pT()/diphoton.Pt());

if (pt_W_reco > 75 /*&& neutrino.pT() >= 20 && gamma1.pT() >=10 && gamma2.pT() >=10 && abs(gamma1.eta()) <= 2.5 && abs(gamma2.eta()) <= 2.5 */) {
h_cos_theta_star_cut1->Fill(cos_theta_star);
h_phi_star_cut1->Fill(phi_star);
}

}
/////////////

if (lepton.pT() >= 30 && abs(lepton.eta()) <= 2.5 && neutrino.pT() >= 20 /*&& gamma1.pT() >=10 && gamma2.pT() >=10 && abs(gamma1.eta()) <= 2.5 && abs(gamma2.eta()) <= 2.5 */) {
h_cos_theta_star_cut2->Fill(cos_theta_star);
h_phi_star_cut2->Fill(phi_star);

}


if (lepton.pT() >= 30 && abs(lepton.eta()) <= 2.5 && gamma1.pT() >=10 && gamma2.pT() >=10 && abs(gamma1.eta()) <= 2.5 && abs(gamma2.eta()) <= 2.5) {
h_cos_theta_star_cut3->Fill(cos_theta_star);
h_phi_star_cut3->Fill(phi_star);
}


if (lepton.pT() >= 30 && abs(lepton.eta()) <= 2.5 && neutrino.pT() >= 20 && gamma1.pT() >=10 && gamma2.pT() >=10 && abs(gamma1.eta()) <= 2.5 && abs(gamma2.eta()) <= 2.5) {
h_cos_theta_star_cut4->Fill(cos_theta_star);
h_phi_star_cut4->Fill(phi_star);
}


tree->Fill();
}} // All Event loop
     
/////////////////////////////////////////////////////////////////////////
//Histogram for cut flow 
    int cut_flow[6] = {n1,n2,n3,n4,n5,n6};
    const char *cuts[6] = {"nEvent","par ev","NJets >=2",">=2 jets with no ALP Gamma","VBF selection","Central Jet Veto"};
    for (int i = 0 ; i < hist_cut_flow->GetNbinsX(); i++){
    hist_cut_flow->SetBinContent(i+1, cut_flow[i]);
    hist_cut_flow->GetXaxis()->SetBinLabel(i+1, cuts[i]);
    }

//////////A_FB////////////////////

TH1F* h_AFB = new TH1F("AFB_ptW", ";p_{T}^{W} [GeV];A_{FB}^{hel}", nAFBbins, wptBins);

for (int i = 0; i < nAFBbins; ++i) {

    double AFB = 0.0;

    if ((Nf[i] + Nb[i]) > 0) {
        AFB = (Nf[i] - Nb[i]) / (Nf[i] + Nb[i]);
    }

    h_AFB->SetBinContent(i+1, AFB);
}


TGraph* g = new TGraph(nAFBbins);

for (int i = 0; i < nAFBbins; ++i) {

    double x = 0.5 * (wptBins[i] + wptBins[i+1]); // bin center
    double y = (Nf[i] + Nb[i] > 0) ? (Nf[i] - Nb[i])/(Nf[i] + Nb[i]) : 0.0;

    g->SetPoint(i, x, y);
}
/////////////////////////////////////////////////////
outfile->Write();
TFile *fout = TFile::Open("/data/SOM_Test/ALP_test/Result_ALP_Wp_125.root","RECREATE");
//TFile *fout = TFile::Open("/data/SOM_Test/ALP_test/Result_Higgs_Wp_125.root","RECREATE");

fout->cd();
hist_cut_flow->Write();
missE_all->Write();
hist_W_pt->Write();
hist_W_pt1->Write();
hist_W_pt2->Write();
hist_W_pt3->Write();
hist_dphi_Wrecopt_digam->Write();
hist_deta_Wrecopt_digam->Write();
hist_dphi_Wrecopt_W->Write();
h2_ptlep_vs_digam->Write();
h2_ptlep_vs_Wreco->Write();
hist_W_eta->Write();
hist_W_mass->Write();
hist_W_phi->Write();
hist_lep_pt->Write();
hist_neutrino_pt->Write();
hist_dphi_lep_MET->Write();
h_cos_theta_star->Write();
h_cos_theta_star_cut1->Write();
h_cos_theta_star_cut2->Write();
h_cos_theta_star_cut3->Write();
h_cos_theta_star_cut4->Write();
h_phi_star_cut->Write();
h_phi_star_cut1->Write();
h_phi_star_cut2->Write();
h_phi_star_cut3->Write();
h_phi_star_cut4->Write();
hist_LP->Write();
h2_cos_vs_ptW->Write();
hist_sig_lep_rap->Write();
hist_dphi_lep_gamma2->Write();
hist_dphi_lep_gamma1->Write();
hist_dR_lep_gamma1->Write();


hist_ratio_ptl_ptW->Write();
hist_ratio_ptl_ptdigam->Write();

h_cos_theta1->Write();

h2_cos_vs_ptlep->Write();
h_AFB->Write();
g->Write("AFB_vs_WpT");

for (int i = 0; i < nPtBins; ++i) {
  h_costheta_wpt[i]->Write();
}

for (int i = 0; i < netaBins; ++i) {
  h_costheta_leta[i]->Write();
}


for (int i = 0; i < nPtBins; ++i) {
  h_costheta_lpt[i]->Write();
}


delete fout;
pythia.stat();
return 0;
}
