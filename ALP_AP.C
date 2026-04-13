#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(/home/shreecheta/project_sc/tools/fastjet-install/lib/libfastjetplugins.so)
//R__LOAD_LIBRARY(libfastjet)
R__LOAD_LIBRARY(/home/shreecheta/project_sc/tools/fastjet-install/lib/libDynamicRPlugin.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "modules/TrackCovariance.h"
#include "external/TrackCovariance/TrkUtil.h"
#include "external/TrackCovariance/VertexFit.h"
#include "external/fastjet/PseudoJet.hh"
#include "external/fastjet/tools/MassDropTagger.hh"
#include "external/fastjet/Selector.hh"
#include "external/fastjet/JetDefinition.hh"
#include "external/fastjet/ClusterSequence.hh"
#include "external/fastjet/contribs/RecursiveTools/SoftDrop.hh" // In external code, this should be fastjet/contrib/SoftDrop.hh
#include "external/fastjet/tools/Filter.hh"
#include "external/fastjet/tools/Pruner.hh"
#include "external/fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "external/fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "external/fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "/home/shreecheta/project_sc/tools/fastjet-3.4.0/DynamicRJetAlgorithm/fastjet/DynamicRJetPlugin.hh"
//------------------------------------------------------------------------------
#else
class ExRootTreeReader;
class ExRootResult;
#endif
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>

using namespace std;
using namespace fastjet;
using namespace contrib;
///////////////////////////////////////////////////////
TLorentzVector reconstructNeutrino(
    const TLorentzVector& lep,
    double met, double met_phi){
    const double MW = 80.379;

    double met_px = met * cos(met_phi);
    double met_py = met * sin(met_phi);

    // --- tunable parameters ---
    double sigma_MET = 10.0;  // Met resolution in GeV 
    double sigma_W   = 3.0;   // W width is 2.1 GeV
    double sigma_pz   = 100.0;   // pz regularization
    //double sigma_pz = 0.5 * (lep.Pt() + met);

    double bestChi2 = 1e20;
    TLorentzVector bestNu;

    // --- scan region (keep this modest) ---
    for (double dx = -10; dx <= 10; dx += 1.0) {
        for (double dy = -10; dy <= 10; dy += 1.0) {

            double nu_px = met_px + dx;
            double nu_py = met_py + dy;

            // --- solve quadratic for pz ---
            double mu = (MW*MW - lep.M2())/2.0
                      + lep.Px()*nu_px + lep.Py()*nu_py;

            double El = lep.E();
            double a  = El*El - lep.Pz()*lep.Pz();
            double b  = -2.0 * mu * lep.Pz();
            double c  = mu*mu - El*El*(nu_px*nu_px + nu_py*nu_py);

            double disc = b*b - 4*a*c;

            if (disc < 0) continue;

            double sqrtD = sqrt(disc);

            for (int sol = 0; sol < 2; sol++) {

                double pz = (sol == 0)
                          ? (-b + sqrtD)/(2*a)
                          : (-b - sqrtD)/(2*a);

                // --- Reject unphysical large pz ---
                //if (fabs(pz) > 300.0) continue;
                          
                double nu_E = sqrt(nu_px*nu_px + nu_py*nu_py + pz*pz);

                TLorentzVector nu(nu_px, nu_py, pz, nu_E);
                TLorentzVector W = lep + nu;

                double mW = W.M();

                // --- χ² ---
                double chi2_MET = (dx*dx + dy*dy)/(sigma_MET*sigma_MET);
                double chi2_W   = (mW - MW)*(mW - MW)/(sigma_W*sigma_W);
                double chi2_pz  = (pz*pz)/(sigma_pz*sigma_pz);

                
                double chi2 = chi2_MET + chi2_W + chi2_pz;
                //double chi2 = chi2_W;

                if (chi2 < bestChi2) {
                    bestChi2 = chi2;
                    bestNu = nu;
                }

                else if (fabs(chi2 - bestChi2) < 0.1) {
                    // tie-breaker: prefer smaller |pz|
                    if (fabs(pz) < fabs(bestNu.Pz())) {
                        bestNu = nu;
                    }
                }    
            }
        }
    }

    // --- fallback if nothing found ---
    if (bestChi2 > 1e10) {
        double mu = (MW*MW - lep.M2())/2.0
                  + lep.Px()*met_px + lep.Py()*met_py;

        double a = lep.Pt()*lep.Pt();
        double b = -2.0 * mu * lep.Pz();

        double pz = -b / (2.0 * a);

        double E = sqrt(met_px*met_px + met_py*met_py + pz*pz);
        bestNu.SetPxPyPzE(met_px, met_py, pz, E);
    }

    return bestNu;
}

////////////////////////////////


void AnalyseEvents(ExRootTreeReader *treeReader)
{

    
  TH1F *hist_lep_rap_gen = new TH1F("lep_rap_gen","lep_rap_gen", 20, -4.5 , 4.5);
  TH1F *hist_mag_3momenta_Wgen = new TH1F("w_3momenta_gen","w_3momenta_gen", 20, 0., 1000);
  TH1F *hist_mag_3momenta_W = new TH1F("w_3momenta","w_3momenta", 20, 0., 1000);
  TH1F *hist_W_mass_gen = new TH1F("W_mass_gen","W_mass_gen", 20, 40., 160);



  TH1F* hist_MET = new TH1F("missing Et", "Missing Et", 100,0.0,1000);
  TH1F* hist_nu_pt = new TH1F("nu_pt", "nu_pt", 100,0.0,1000);
  TH1F* hist_nu_pz_diff_gen = new TH1F("nu_pz_diff_gen", "nu_pz_gen", 15,-10,140);
  TH1F* hist_nu_pz_diff = new TH1F("nu_pz_diff", "nu_pz_diff", 15,-10,140);

  TH1F* hist_nu_pxpy_diff_gen = new TH1F("nu_pxpy_diff_gen", "nu_pxpy_diff_gen", 15,-10,140);
  TH1F* hist_nu_pxpy_diff = new TH1F("nu_pxpy_diff", "nu_pxpy_diff", 15,-10,140);


  TH1F* hist_MET_gen = new TH1F("missing_Et_gen", "Missing_Et_gen", 100,0.0,1000);
  TH1F *hist_dphi_Wrecopt_digam_gen = new TH1F("dphi_Wrecopt_digam_gen","dphi_Wrecopt_digam_gen", 100, -6.0, 6.);

  TH1F* hist_DR_gen_Del_lep = new TH1F("DR_gen_Del_lep", "DR_gen_Del_lep", 100,0.0,1);

  TH1F *hist_dphi_lep_MET_gen = new TH1F("dphi_lep_MET_gen","dphi_lep_MET_gen", 100, -6.0, 6.);

  TH1F* hist_neutrino_pt = new TH1F("neutrino_pt", "neutrino_pt", 100,0.0,1000);

  TH1F* hist_Ngam = new TH1F("Ngam", "Ngam", 10,0.0,10);
  TH1F* hist_lep_size = new TH1F("Nlep", "Nlep",5,0,5);

  TH1F* hist_inv_diphoton = new TH1F("inv_diphoton", "inv_diphoton", 100,0.0,100);
  TH1F *hist_ratio_ptl_ptW = new TH1F("ratio_ptl_ptW","ratio_ptl_ptW", 100, 0, 1);
  TH1F *hist_ratio_ptl_ptdigam = new TH1F("ratio_ptl_ptdigam","ratio_ptl_ptdigam", 100, 0, 1);

  //TH1F *hist_= new TH1F("delta R of two gamma from ALP","delta R of two gamma from ALP", 40, 0.0, 1.);
  TH1F *hist_W_pt_gen = new TH1F("pt_W_gen","pt_W_gen", 100, 0.0, 1000.);
  TH1F *hist_W_pt1 = new TH1F("pt_W1","pt_W1", 100, 0.0, 1000.);
  TH1F *hist_W_pt2 = new TH1F("pt_W2","pt_W2", 100, 0.0, 1000.);
  TH1F *hist_W_pt3 = new TH1F("pt_W3","pt_W3", 100, 0.0, 1000.);
  TH1F *hist_dphi_Wrecopt_digam = new TH1F("dphi_Wrecopt_digam","dphi_Wrecopt_digam", 100, -6.0, 6.);
  TH1F *hist_deta_Wrecopt_digam = new TH1F("deta_Wrecopt_digam","deta_Wrecopt_digam", 100, -6, 6.);

  TH1F *hist_dphi_Wrecopt_W = new TH1F("dphi_Wrecopt_W","dphi_Wrecopt_W", 100, -6.0, 6.);
  TH2F* h2_ptlep_vs_Wreco = new TH2F("ptlep_vs_Wreco", "p_{T}^{l} vs reco p_{T}^{W}; reco p_{T}^{W} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  TH2F* h2_ptlep_vs_digam = new TH2F("ptlep_vs_digam", "p_{T}^{l} vs p_{T}^{2#gamma}; p_{T}^{2#gamma} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  TH1F *hist_LP = new TH1F("LP","LP", 100, -2, 2.);
  TH1F *hist_LP_gen = new TH1F("LP_gen","LP_gen", 100, -2, 2.);
  

  TH1F *hist_W_eta = new TH1F("eta_W","eta_W", 100, -6.0, 6.);
  TH1F *hist_W_mass = new TH1F("W_mass","W_mass", 20, 40., 160);

  TH1F *hist_W_phi = new TH1F("phi_W","phi_W", 100, -6.0, 6.);
  TH1F* h_cos_theta_star = new TH1F("cos_theta_star", "cos_theta_star", 20, -1.0, 1.0); 

  TH1F* h_cos_theta_star_from_LP = new TH1F("cos_theta_star_from_LP", "cos_theta_star_from_LP", 20, -1.0, 1.0);  

  TH1F* h_cos_theta_star_cut1_gen = new TH1F("cos_theta_star_cut1_gen", "cos_theta_star_cut1_gen", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_cut1 = new TH1F("cos_theta_star_cut1", "cos_theta_star_cut1", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_from_LP_cut1 = new TH1F("cos_theta_star_from_LP_cut1", "cos_theta_star_from_LP_cut1", 20, -1.0, 1.0);  
  TH1F *hist_LP_cut1 = new TH1F("LP_cut1","LP_cut1", 100, -2, 2.);
  TH1F *hist_W_pt3_cut1 = new TH1F("pt_W3_cut1","pt_W3_cut1", 100, 0.0, 1000.);
  TH1F *hist_lep_pt_cut1 = new TH1F("lep_pt_cut1","lep_pt_cut1", 25, 0.0, 250.);
  TH1F *hist_dphi_lep_MET_cut1 = new TH1F("dphi_lep_MET_cut1","dphi_lep_MET_cut1", 100, -6 , 6);
  TH1F *hist_ratio_ptl_ptW_cut1 = new TH1F("ratio_ptl_ptW_cut1","ratio_ptl_ptW_cut1", 100, 0, 1);

  TH1F* h_cos_theta_star_gen = new TH1F("cos_theta_star_gen", "cos_theta_star_gen", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_new_gen = new TH1F("cos_theta_star_new_gen", "cos_theta_star_new_gen", 20, 0., 1.0); 
  TH1F* h_cos_theta_star_new = new TH1F("cos_theta_star_new", "cos_theta_star_new", 20, 0., 1.0); 

  TH1F* h_cos_theta_star_new_cut1 = new TH1F("cos_theta_star_new_cut1", "cos_theta_star_new_cut1", 20, 0., 1.0); 
  TH1F* h_cos_theta_star_new_cut1_gen = new TH1F("cos_theta_star_new_cut1_gen", "cos_theta_star_new_cut1_gen", 20, 0., 1.0); 


  
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
  TH1F *hist_lep_pt_gen = new TH1F("lep_pt_gen","lep_pt_gen", 25, 0.0, 250.);

  TH1F *hist_dphi_lep_MET = new TH1F("dphi_lep_MET","dphi_lep_MET", 100, -6 , 6);
  TH1F *hist_sig_lep_rap = new TH1F("sig_lep_rap","sig_lep_rap", 20, -4.5 , 4.5);


  /////////////////////////////////////////
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  //TClonesArray *branchEflow = treeReader->UseBranch("Eflow");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack"); // Changed from Eflow to EFlowTrack
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron"); // Added EFlowNeutralHadron
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton"); // Added EFlowPhoton
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGENMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchTower = treeReader->UseBranch("Tower"); // Added Tower branch explicitly
  
  
  //=========Declare=========================
  Long64_t allEntries = treeReader->GetEntries();
  Double_t m_h = 0;
  cout << "Chain contains " << allEntries << " events" << endl;

  Electron *electron;
  Photon *photon;
  Muon *muon;
  MissingET *met;
  MissingET *Genmet;
  ScalarHT *HT;
  Track *track;
  Tower *tower;
  Jet *jet;
  //TObject *object;
  //TLorentzVector jetMomentum;
  GenParticle *particle;
  GenParticle *mother;
  Int_t i, j;
  //////////////////Varibles ///////
   //TFile *outfile = new TFile("/data/SOM_Test/ALP_test/OutputFile_Higgs_Wp_125_Delphes.root","RECREATE");
   TFile *outfile = new TFile("/data/SOM_Test/ALP_test/OutputFile_Higgs_Wp_125_Delphes.root","RECREATE");

  TTree *tree=new TTree("Tree","Signal");
  Float_t  W_reco_pt = 0., Lep_pt=0., MET_et = 0.0,  Cos_theta_star =0.0, Dphi_Lep_MET=0., Dphi_Wrecp_Diphoton=0.,Deta_Wrecp_Diphoton=0., Phi_star=0.;
  Float_t LP=0.;
  Int_t Total_Events=0 ;

  tree->Branch("Total_Events" ,     & Total_Events ,      "Total_Events/I");
  tree->Branch("W_reco_pt",&W_reco_pt,"W_reco_pt/F");
  tree->Branch("Lep_pt",&Lep_pt,"Lep_pt/F");
  tree->Branch("MET_et",&MET_et,"MET_et/F");
  tree->Branch("Dphi_Lep_MET",&Dphi_Lep_MET,"Dphi_Lep_MET/F");
  tree->Branch("Cos_theta_star",&Cos_theta_star,"Cos_theta_star/F");
  tree->Branch("Deta_Wrecp_Diphoton",&Deta_Wrecp_Diphoton,"Deta_Wrecp_Diphoton/F");
  tree->Branch("Phi_star",&Phi_star,"Phi_star/F");
  tree->Branch("Dphi_Wrecp_Diphoton",&Dphi_Wrecp_Diphoton,"Dphi_Wrecp_Diphoton/F");
  tree->Branch("LP",&LP,"LP/F");

//=================================
  Long64_t entry;

  int n1 = 0 , n2 = 0 , n3 = 0, n4 =0 , n5 = 0 , n6 =0, n7=0, n8=0,n9=0,n10=0,n11=0,n12=0,n13=0, n14=0, n15=0;   
  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry) //allEntries
  { 
     int x_entry = allEntries/10;
    //if (entry < 5700) continue;
    //if (entry > 2) continue;
    //if (entry  != 53) continue;

    if(entry%x_entry == 0) cout << "event number: " << entry << endl;
    
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    // test Particle branch
    Int_t Ngen = branchParticle->GetEntries();
    Total_Events = allEntries;
//======================Preselection loop========================== 

  //=======================Analyse missing ET ============================================================================
if(branchMissingET->GetEntriesFast() > 0){
      met = (MissingET*) branchMissingET->At(0);
} 
TVector2 Pmet;
Pmet.SetMagPhi(met->MET, met->Phi);

if(branchGENMissingET->GetEntriesFast() > 0){
      Genmet = (MissingET*) branchGENMissingET->At(0);
} 

//hist_MET->Fill(met->MET);

TLorentzVector MISSE = met->P4();
//cout<<MISSE.Pt()<<"    "<<met->MET<<"   "<<MISSE.Pz()<<"    "<<MISSE.E()<<endl;
// cout<<TMath::Sqrt(MISSE.Px() * MISSE.Px()  +      MISSE.Py() * MISSE.Py() + MISSE.Pz() * MISSE.Pz() ) <<endl;
//cout<<MISSE.Px()<<"   "<<MISSE.Py()<<"   "<<MISSE.Pz()<<"   "<<MISSE.E()<<endl;

fastjet::PseudoJet METlike (MISSE.Px(), MISSE.Py(), MISSE.Pz(), MISSE.E());

int W_idx = -1;
int lepton_idx = -1;
int neutrino_idx = -1;
int ALP_idx = -1;
int light_dau1 = -1;
int light_dau2 = -1;
int ALP_gam1 = -1;
int ALP_gam2 = -1;
////////////////////////////////////////////////////////////////////
for(int i = 0; i < branchParticle->GetEntriesFast(); ++i) {
    particle = (GenParticle*) branchParticle->At(i);

    if (abs (particle->PID) == 24 && abs(particle->Status) == 62){
        W_idx = i;
        //cout<<particle->Status<<"  "<<entry<<endl;
    }
    if ( (abs(particle->PID) == 25 ||  abs(particle->PID) == 9000005 ||  abs(particle->PID) == 36) && abs(particle->Status) == 62) {  // ALP
                ALP_idx= i;
    }
    if (abs(particle->PID) == 24 && abs(particle->Status) == 62) {  // W boson (after FSR)
                W_idx = i;
                // Check daughters
                int d1 = particle->D1;
                int d2 = particle->D2;
                
                if (d1 > 0 && d2 > 0) {

                    auto* daughter1= (GenParticle*) branchParticle->At(d1);
                    auto* daughter2= (GenParticle*) branchParticle->At(d2);
                    
                    // Check for leptonic decay (e, mu, tau)
                    if (abs(daughter1->PID) == 11 || abs(daughter1->PID) == 13 || abs(daughter1->PID) == 15) {
                        lepton_idx = d1;
                        neutrino_idx = d2;

                          //cout<<d1<<"  "<<d2<<" "<<iEvent<<endl;
                   
                    } else if (abs(daughter2->PID) == 11 || abs(daughter2->PID) == 13 || abs(daughter2->PID) == 15) {
                        lepton_idx = d2;
                        neutrino_idx = d1;
                    
                    }

                }
    }
    if (abs(particle->PID) == 25 || abs(particle->PID) == 9000005) {
  
    auto* daughter1= (GenParticle*) branchParticle->At(particle->D1);
    auto* daughter2= (GenParticle*) branchParticle->At(particle->D2);

    bool daughters_valid = true;
    
    if (abs(particle->PID) == 25 && (abs(daughter1->PID) == 25 || abs(daughter2->PID) == 25)) {
        daughters_valid = false;
    }
    if (abs(particle->PID) == 9000005 && (abs(daughter1->PID) == 9000005 || abs(daughter2->PID) == 9000005)) {
        daughters_valid = false;
    }

    if(daughters_valid == false) continue;

   if (abs(daughter1->PID) == 22) {
            ALP_gam1 = particle->D1;
        }

   if (abs(daughter2->PID) == 22){ 
            ALP_gam2 = particle->D2;
      }
    }
} //end of genparticle loop
////////////////////////Gen Level Checks/////////////////////
GenParticle *W = (GenParticle*) branchParticle->At(W_idx);
GenParticle *ALP = (GenParticle*) branchParticle->At(ALP_idx);
GenParticle *lepton = (GenParticle*) branchParticle->At(lepton_idx); 
GenParticle *neutrino = (GenParticle*) branchParticle->At(neutrino_idx); 
GenParticle *gamma1 = (GenParticle*) branchParticle->At(ALP_gam1); 
GenParticle *gamma2 = (GenParticle*) branchParticle->At(ALP_gam2);
TLorentzVector lep_vec_org_gen;
if (lepton->PT > 20 && abs(lepton->Eta) < 2.5  && gamma1->PT > 10 && gamma2->PT > 10 && abs(gamma1->Eta) < 2.5 && 
abs(gamma2->Eta) < 2.5  && Genmet->MET > 20){ 

    TLorentzVector W_vec_gen(W->Px, W->Py, W->Pz, W->E);
    TLorentzVector neutrino_vec_gen(neutrino->Px, neutrino->Py, neutrino->Pz, neutrino->E);
    TLorentzVector lep_vec_gen(lepton->Px, lepton->Py, lepton->Pz, lepton->E);
    TLorentzVector gam1_vec_gen(gamma1->Px, gamma1->Py, gamma1->Pz, gamma1->E);
    TLorentzVector gam2_vec_gen(gamma2->Px, gamma2->Py, gamma2->Pz, gamma2->E);
    TLorentzVector diphoton_gen = gam1_vec_gen + gam2_vec_gen ;
    lep_vec_org_gen.SetPxPyPzE(0,0,0,0);
    lep_vec_org_gen = lep_vec_gen;
    TLorentzVector nu_reco_gen(0,0,0,0);
    nu_reco_gen = reconstructNeutrino(lep_vec_org_gen, Genmet->MET, Genmet->Phi);
    //cout<<nu_reco.Pt()<<"   "<<W_reco.M()<<endl;
    // Reconstruct W
    TLorentzVector W_reco_gen = lep_vec_org_gen + nu_reco_gen;

    if(nu_reco_gen.Pt() == 0 ) continue;

    //cout<<W_vec_gen.Vect().Mag()<<"  "<<W_reco_gen.Vect().Mag()<<endl;
    hist_W_mass_gen->Fill(W_reco_gen.M());
    hist_MET_gen->Fill(Genmet->MET);
    hist_lep_pt_gen->Fill(lepton->PT);
    hist_lep_rap_gen->Fill(lep_vec_org_gen.Eta());
    hist_nu_pt->Fill(neutrino->PT);

    hist_nu_pz_diff_gen->Fill(abs(neutrino->Pz) - abs(nu_reco_gen.Pz()));
    hist_nu_pxpy_diff_gen->Fill(sqrt(pow(neutrino->Px - nu_reco_gen.Px(),2) + sqrt(pow(neutrino->Py - nu_reco_gen.Py(),2))));


    hist_dphi_lep_MET_gen->Fill(abs(lep_vec_org_gen.DeltaPhi(neutrino_vec_gen)));
    hist_dphi_Wrecopt_digam_gen->Fill(abs(W_vec_gen.DeltaPhi(diphoton_gen)));

    hist_mag_3momenta_Wgen->Fill(W_vec_gen.Vect().Mag());
    hist_W_pt_gen->Fill(W_vec_gen.Pt());
    
    TVector3 boostW = W_vec_gen.BoostVector();
    lep_vec_gen.Boost(-boostW);
    TVector3 lep_dir_Wrest = lep_vec_gen.Vect().Unit();   // lepton in W rest frame
    TVector3 W_dir_lab     = W_vec_gen.Vect().Unit();     // W direction in lab
    double cos_theta_star_gen = lep_dir_Wrest.Dot(W_dir_lab);
    h_cos_theta_star_gen->Fill(cos_theta_star_gen);
    double costhetastar_equivalant = abs(lepton->E - neutrino->E)/W_vec_gen.Vect().Mag();

    h_cos_theta_star_new_gen->Fill(costhetastar_equivalant);
    
    /////LP 
    double LP_gen= (lep_vec_org_gen.Px()*W_vec_gen.Px() + lep_vec_org_gen.Py()*W_vec_gen.Py())
                  / (W_vec_gen.Px()*W_vec_gen.Px() + W_vec_gen.Py()*W_vec_gen.Py());

    hist_LP_gen->Fill(LP_gen);

    if (W_vec_gen.Pt() > 75) {
        h_cos_theta_star_cut1_gen->Fill(cos_theta_star_gen);
        h_cos_theta_star_new_cut1_gen->Fill(costhetastar_equivalant);
        
    }

}
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//vector<Photon*> Photon_vec;
std::vector<fastjet::PseudoJet> Photon_vec;

if (branchPhoton->GetEntriesFast() > 1){
for(int k = 0; k < branchPhoton->GetEntriesFast(); ++k ){
 photon = (Photon*) branchPhoton->At(k);     
 TLorentzVector photonMomentum = photon->P4();
 if(photonMomentum.Pt() > 10 && fabs(photonMomentum.Eta())< 2.5) {   
    fastjet::PseudoJet pseudoJet(photonMomentum.Px(), photonMomentum.Py(), 
                                     photonMomentum.Pz(), photonMomentum.E());
    Photon_vec.push_back(pseudoJet);

        }
    } 
    hist_Ngam->Fill(Photon_vec.size());
    //cout<<Photon_vec.size()<<endl;
}
if (Photon_vec.size()> 2){
    Photon_vec = sorted_by_pt(Photon_vec);
    hist_inv_diphoton->Fill( (Photon_vec[0] + Photon_vec[1]).m() );
  }


std::vector<fastjet::PseudoJet> Lepton_vec;

for(Int_t i = 0; i < branchElectron->GetEntriesFast(); ++i) {
    Electron *electron = (Electron*) branchElectron->At(i);
    if(electron->PT > 20.0 && fabs(electron->Eta) < 2.5 ) {
         TLorentzVector leptonMomentum = electron->P4();
        fastjet::PseudoJet pseudoJet(leptonMomentum.Px(), leptonMomentum.Py(), 
                                     leptonMomentum.Pz(), leptonMomentum.E());
        Lepton_vec.push_back(pseudoJet);
      }
    }
for(Int_t i = 0; i < branchMuon->GetEntriesFast(); ++i) {
    Muon *muon = (Muon*) branchMuon->At(i);
    if(muon->PT > 20.0 && fabs(muon->Eta) < 2.5) {
        TLorentzVector leptonMomentum = muon->P4();
        fastjet::PseudoJet pseudoJet(leptonMomentum.Px(), leptonMomentum.Py(), 
                                     leptonMomentum.Pz(), leptonMomentum.E());
        Lepton_vec.push_back(pseudoJet);
      }
    }
//cout<<Lepton_vec.size()<<endl;
hist_lep_size->Fill(Lepton_vec.size());

////////////////////////////////////Delphes Level///////////////////////////
double cos_theta_star =-5, phi_star = -5;    
if(Lepton_vec.size() >= 1 && Photon_vec.size() >= 2 && met->MET > 20){ //////////////////////////Basic Selection

    Lepton_vec = sorted_by_pt(Lepton_vec);
    //cout<<Photon_vec[0].delta_R(Photon_vec[1])<<"  "<<Lepton_vec[0].delta_R(Photon_vec[1])<<endl;

    TLorentzVector lep_vec(Lepton_vec[0].px(), Lepton_vec[0].py(), Lepton_vec[0].pz(), Lepton_vec[0].e());
    TLorentzVector lep_vec_org = lep_vec;
    // Reconstruct neutrino
    TLorentzVector nu_reco = reconstructNeutrino(lep_vec_org, met->MET, met->Phi);
    TLorentzVector W_reco = lep_vec_org + nu_reco;

    //cout<<nu_reco.Pt()<<"   "<<W_reco.M()<<endl;
    if(nu_reco.Pt() == 0 ) continue;

    hist_nu_pz_diff->Fill(abs(neutrino->Pz) - abs(nu_reco.Pz()));
    hist_nu_pxpy_diff->Fill(sqrt(pow(neutrino->Px - nu_reco.Px(),2) + sqrt(pow(neutrino->Py - nu_reco.Py(),2))));

    hist_W_mass->Fill(W_reco.M());


    hist_MET->Fill(met->MET);
    hist_DR_gen_Del_lep->Fill(lep_vec_org.DeltaR(lep_vec_org_gen));

    TLorentzVector gam1_vec(Photon_vec[0].px(), Photon_vec[0].py(), Photon_vec[0].pz(), Photon_vec[0].e());
    TLorentzVector gam2_vec(Photon_vec[1].px(), Photon_vec[1].py(), Photon_vec[1].pz(), Photon_vec[1].e());
    TLorentzVector diphoton = gam1_vec + gam2_vec;
    hist_W_pt1->Fill(diphoton.Pt());
    
    // Reconstruct W
    hist_dphi_lep_MET->Fill(abs(lep_vec_org.DeltaPhi(nu_reco)));

    hist_mag_3momenta_W->Fill(W_reco.Vect().Mag());

    //cout<<W->PT<<"   "<<W_reco.Pt()<<"   "<<W_reco.M()<<endl;

    W_reco_pt = 0., Lep_pt = 0., MET_et=0., Dphi_Wrecp_Diphoton=0., Deta_Wrecp_Diphoton=0.;
    W_reco_pt = W_reco.Pt();
    Lep_pt = lep_vec_org.Pt();
    MET_et = met->MET;
    Dphi_Wrecp_Diphoton=abs(W_reco.DeltaPhi(diphoton));
    Deta_Wrecp_Diphoton=abs(W_reco.Eta() - diphoton.Eta()); 

    
    hist_W_pt2->Fill((MISSE + lep_vec_org).Pt()); // Vec sum of MET & Lep 
    hist_W_pt3->Fill(W_reco.Pt()); /////// var 
    
    hist_dphi_Wrecopt_digam->Fill(Dphi_Wrecp_Diphoton);
    hist_deta_Wrecopt_digam->Fill(Deta_Wrecp_Diphoton);
    //hist_dphi_Wrecopt_W->Fill(W_reco.DeltaPhi(W_vec));
    h2_ptlep_vs_Wreco->Fill(lep_vec_org.Pt(), W_reco.Pt());
    h2_ptlep_vs_digam->Fill(lep_vec_org.Pt(), diphoton.Pt());

    /////////Costheta*   
    //TLorentzVector W_vec(W->Px, W->Py, W->Pz, W->E);  
    //TVector3 boostW = W_vec.BoostVector();
    TVector3 boostW = W_reco.BoostVector();
    lep_vec.Boost(-boostW);

    TVector3 lep_dir_Wrest = lep_vec.Vect().Unit();   // lepton in W rest frame
    TVector3 W_dir_lab     = W_reco.Vect().Unit();     // W direction in lab

    cos_theta_star = lep_dir_Wrest.Dot(W_dir_lab);

    Cos_theta_star =0, Phi_star = 0;
    Cos_theta_star = cos_theta_star;
    h_cos_theta_star->Fill(cos_theta_star);    

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
    h_cos_theta1->Fill(cos_theta1);

    // check 
    //cout<<" In CM this should be -1 => "<<W_dir_CM.Dot(H_dir_CM)<<endl;

//////Phi*
//phi_star azimuthal angle of the charged lepton in the W rest frame
//measured with respect to the production plane
TVector3 beam_dir(0.0, 0.0, 1.0);   // proton beam direction
TVector3 y_hat = beam_dir.Cross(W_dir_lab).Unit();
TVector3 x_hat = y_hat.Cross(W_dir_lab).Unit();

phi_star = abs(atan2(lep_dir_Wrest.Dot(y_hat),lep_dir_Wrest.Dot(x_hat)));
Phi_star = phi_star;
h_phi_star_cut->Fill(Phi_star);

hist_lep_pt->Fill(lep_vec_org.Pt());
hist_sig_lep_rap->Fill(lep_vec_org.Eta());

hist_dphi_lep_gamma1->Fill(lep_vec_org.DeltaPhi(gam1_vec));
hist_dphi_lep_gamma2->Fill(lep_vec_org.DeltaPhi(gam2_vec));
hist_dR_lep_gamma1->Fill(lep_vec_org.DeltaR(gam2_vec));

Dphi_Lep_MET=0.;
Dphi_Lep_MET=abs(lep_vec_org.DeltaPhi(nu_reco));

//////LP
LP=0.;
TVector3 pT_3lep = lep_vec_org.Vect();
TVector3 pT_3W   = W_reco.Vect();

//LP= pT_3lep.Dot(pT_3W)/pT_3W.Mag2();
LP = (lep_vec_org.Px()*W_reco.Px() + lep_vec_org.Py()*W_reco.Py())
                  / (W_reco.Px()*W_reco.Px() + W_reco.Py()*W_reco.Py());

h_cos_theta_star_new->Fill(abs(lep_vec_org.E() - nu_reco.E())/pT_3W.Mag());

//cout<<cos_theta_star<<"   "<<2*(LP - 0.5)<<"  "<< abs(lep_vec_org.E() - MET_et)/W_reco.Pt()<<endl;

hist_LP->Fill(LP);
h_cos_theta_star_from_LP->Fill(2*(LP - 0.5));
//if (W_reco.Pt() > 50) cout<<Cos_theta_star<<"   "<<2*(LP - 0.5)<<"   "<<entry<<endl;

hist_ratio_ptl_ptW->Fill(lep_vec_org.Pt()/W_reco.Pt());
hist_ratio_ptl_ptdigam->Fill(lep_vec_org.Pt()/diphoton.Pt());

if(W_reco.Pt() > 75){
    //cout<<W_reco.M()<<endl;
    h_cos_theta_star_from_LP_cut1->Fill(2*(LP - 0.5));
    h_cos_theta_star_cut1->Fill(Cos_theta_star);
    hist_LP_cut1->Fill(LP);
    hist_W_pt3_cut1->Fill(W_reco.Pt()); /////// var 
    hist_lep_pt_cut1->Fill(lep_vec_org.Pt());
    hist_ratio_ptl_ptW_cut1->Fill(lep_vec_org.Pt()/W_reco.Pt());
    hist_dphi_lep_MET_cut1->Fill(abs(lep_vec_org.DeltaPhi(nu_reco)));
    h_cos_theta_star_new_cut1->Fill(abs(lep_vec_org.E() - nu_reco.E())/pT_3W.Mag());

}
}








tree->Fill();

  } //All event loop
cout<<n1<<"  "<<n2<<"  "<<n3<<"  "<<n4<<"   "<<n5<<"  "<<n6<<"   "<<n7<<"   "<<n8<<"   "<<n9<<"  "<<n10<<"  "<<
n11<<"   "<<n12<<"  "<<n13<<"   "<<n14<<"  "<<n15<<endl;
//====================================================================
    int cut_flow[14] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14};
    const char *cuts[14] = {"nEvent","par ev","ev", "NJets >=2","VBF selection",">= 2 iso photons","N_FJ >= 1","DR(FJ,VBF Jet) >=2.5","Grooming",
    "tai21 < 0.3","Ch_trk <3", "EFrac gam > 0.7", "N_subjet >=2","lep veto"};
    for (int i = 0 ; i < hist_cut_flow->GetNbinsX(); i++){
    hist_cut_flow->SetBinContent(i+1, cut_flow[i]);
    hist_cut_flow->GetXaxis()->SetBinLabel(i+1, cuts[i]);
    }
outfile->Write();
 //TFile *fout = TFile::Open("/data/SOM_Test/ALP_test/Result_Higgs_Wp_125_Delphes.root","RECREATE");
 TFile *fout = TFile::Open("/data/SOM_Test/ALP_test/Result_Higgs_Wp_125_Delphes.root","RECREATE");

fout->cd();

hist_W_mass_gen->Write();
hist_W_mass->Write();
hist_nu_pt->Write();
hist_nu_pz_diff_gen->Write();
hist_nu_pz_diff->Write();
hist_nu_pxpy_diff_gen->Write();
hist_nu_pxpy_diff->Write();

hist_mag_3momenta_Wgen->Write();
hist_mag_3momenta_W->Write();
hist_DR_gen_Del_lep->Write();
hist_lep_rap_gen->Write();
hist_dphi_Wrecopt_digam_gen->Write();
hist_dphi_lep_MET_gen->Write();
hist_lep_pt_gen->Write();
hist_LP_gen->Write();
hist_MET_gen->Write();
h_cos_theta_star_gen->Write();
h_cos_theta_star_new_gen->Write();
hist_MET->Write();
hist_Ngam->Write();
hist_lep_size->Write();
hist_inv_diphoton->Write();
hist_W_pt_gen->Write();
hist_W_pt1->Write();
hist_W_pt2->Write();
hist_W_pt3->Write();
hist_dphi_Wrecopt_digam->Write();
hist_deta_Wrecopt_digam->Write();
hist_dphi_Wrecopt_W->Write();
h2_ptlep_vs_digam->Write();
h2_ptlep_vs_Wreco->Write();
h_cos_theta_star->Write();
h_cos_theta_star_from_LP->Write();
h_phi_star_cut->Write();
hist_ratio_ptl_ptW->Write();
hist_ratio_ptl_ptdigam->Write();
hist_LP->Write();
h_cos_theta_star_new->Write();
h_cos_theta_star_cut1_gen->Write();
h_cos_theta_star_new_cut1->Write();
h_cos_theta_star_new_cut1_gen->Write();
hist_lep_pt->Write();
hist_sig_lep_rap->Write();

hist_dphi_lep_MET->Write();
hist_dphi_lep_gamma2->Write();
hist_dR_lep_gamma1->Write();

//cut1
h_cos_theta_star_from_LP_cut1->Write();
h_cos_theta_star_cut1->Write();
hist_LP_cut1->Write();
hist_ratio_ptl_ptW_cut1->Write();
hist_W_pt3_cut1->Write();
hist_lep_pt_cut1->Write();
hist_dphi_lep_MET_cut1->Write();
h_cos_theta1->Write();




hist_cut_flow->Write();
delete fout;

}// End of AnalyseEvents

void ALP_AP(const char *inputFile)
{
    gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 // ExRootResult *result = new ExRootResult();
  AnalyseEvents(treeReader);

  cout << "** Exiting..." << endl;

  delete treeReader;
  delete chain;
}
