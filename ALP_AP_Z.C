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
////////////////////////////////



//------------------------------------------------------------------------------


void AnalyseEvents(ExRootTreeReader *treeReader)
{

    
  TH1F *hist_lep_rap_gen = new TH1F("lep_rap_gen","lep_rap_gen", 40, -4.5 , 4.5);
  TH1F *hist_mag_3momenta_Wgen = new TH1F("w_3momenta_gen","w_3momenta_gen", 20, 0., 1000);
  TH1F *hist_mag_3momenta_W = new TH1F("w_3momenta","w_3momenta", 20, 0., 1000);
  TH1F *hist_W_mass_gen = new TH1F("W_mass_gen","W_mass_gen", 200, 40., 160);
  TH1F *hist_W_mass_gen1 = new TH1F("W_mass_gen1","W_mass_gen1", 200, 40., 160);

  TH1F *hist_W_rap = new TH1F("W_rap","W_rap", 40, -4.5 , 4.5);
  TH1F *hist_W_rap_gen = new TH1F("W_rap_gen","W_rap_gen", 40, -4.5 , 4.5);
  TH1F *hist_W_rap_gen1 = new TH1F("W_rap_gen1","W_rap_gen1", 40, -4.5 , 4.5);



  TH1F* hist_MET = new TH1F("missing Et", "Missing Et", 100,0.0,1000);
  TH1F* hist_nu_pt = new TH1F("nu_pt", "nu_pt", 100,0.0,1000);
  TH1F* hist_nu_pz_gen = new TH1F("nu_pz_gen", "nu_pz_gen", 100,0.0,1000);

  TH1F* hist_nu_pz_diff_gen = new TH1F("nu_pz_diff_gen", "nu_pz_gen", 15,-10,140);
  TH1F* hist_nu_pz_diff = new TH1F("nu_pz_diff", "nu_pz_diff", 15,-10,140);

  TH1F* hist_nu_pxpy_diff_gen = new TH1F("nu_pxpy_diff_gen", "nu_pxpy_diff_gen", 15,-10,140);
  TH1F* hist_nu_pxpy_diff = new TH1F("nu_pxpy_diff", "nu_pxpy_diff", 15,-10,140);


  TH1F* hist_MET_gen = new TH1F("missing_Et_gen", "Missing_Et_gen", 100,0.0,1000);
  TH1F *hist_dphi_Wrecopt_digam_gen = new TH1F("dphi_Wrecopt_digam_gen","dphi_Wrecopt_digam_gen", 100, -6.0, 6.);
  TH1F *hist_deta_Wrecopt_digam_gen = new TH1F("deta_Wrecopt_digam_gen","deta_Wrecopt_digam_gen", 100, -6.0, 6.);

  TH1F* hist_DR_gen_Del_lep = new TH1F("DR_gen_Del_lep", "DR_gen_Del_lep", 100,0.0,1);

  TH1F *hist_dphi_lep_MET_gen = new TH1F("dphi_lep_MET_gen","dphi_lep_MET_gen", 100, -6.0, 6.);
  TH1F *hist_deta_lep_MET_gen = new TH1F("deta_lep_MET_gen","deta_lep_MET_gen", 100, -6.0, 6.);


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
  TH2F* h2_ptlep_vs_Wreco = new TH2F("ptlep_vs_Wreco", "p_{T}^{l} vs reco p_{T}^{Z}; reco p_{T}^{Z} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  TH2F* h2_ptlep_vs_digam = new TH2F("ptlep_vs_digam", "p_{T}^{l} vs p_{T}^{2#gamma}; p_{T}^{2#gamma} [GeV]; p_{T}^{l} [GeV]", 100, 0, 1000, 100, 0., 1000);
  TH1F *hist_LP = new TH1F("LP","LP", 50, -2, 2.);
  TH1F *hist_LP_gen = new TH1F("LP_gen","LP_gen", 50, -2, 2.);
  

  TH1F *hist_W_eta = new TH1F("eta_W","eta_W", 100, -6.0, 6.);
  TH1F *hist_W_mass = new TH1F("W_mass","W_mass", 200, 40., 160);

  TH1F *hist_W_phi = new TH1F("phi_W","phi_W", 100, -6.0, 6.);
  TH1F* h_cos_theta_star = new TH1F("cos_theta_star", "cos_theta_star", 20, -1.0, 1.0); 

  TH1F* h_cos_theta_star_from_LP = new TH1F("cos_theta_star_from_LP", "cos_theta_star_from_LP", 20, -1.0, 1.0);  

  TH1F* h_cos_theta_star_cut1_gen = new TH1F("cos_theta_star_cut1_gen", "cos_theta_star_cut1_gen", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_cut1 = new TH1F("cos_theta_star_cut1", "cos_theta_star_cut1", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_from_LP_cut1 = new TH1F("cos_theta_star_from_LP_cut1", "cos_theta_star_from_LP_cut1", 20, -1.0, 1.0);  
  TH1F *hist_LP_cut1 = new TH1F("LP_cut1","LP_cut1", 50, -2, 2.);
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

  TH2F* h2_cos_vs_ptW = new TH2F("cos_vs_ptW", "cos#theta* vs p_{T}^{Z};p_{T}^{Z} [GeV];cos#theta*", 30, 0, 300, 20, -1, 1);
  TH2F* h2_cos_vs_ptlep = new TH2F("cos_vs_ptlep", "cos#theta* vs p_{T}^{lep};p_{T}^{lep} [GeV];cos#theta*", 10, 0, 100, 20, -1, 1);

  TH1F *hist_dphi_lep_gamma1 = new TH1F("dphi_lep_gamma1","dphi_lep_gamma1", 20, 0.0 , 3.5);
  TH1F *hist_dphi_lep_gamma2 = new TH1F("dphi_lep_gamma2","dphi_lep_gamma2", 20, 0.0 , 3.5);
  TH1F *hist_dR_lep_gamma1 = new TH1F("dR_lep_gamma1","dR_lep_gamma1", 20, 0.0 , 3.5);



  TH1F *hist_lep_pt = new TH1F("lep_pt","lep_pt", 25, 0.0, 250.);
  TH1F *hist_lep_pt_gen = new TH1F("lep_pt_gen","lep_pt_gen", 25, 0.0, 250.);

  TH1F *hist_dphi_lep_MET = new TH1F("dphi_lep_MET","dphi_lep_MET", 100, -6 , 6);
  TH1F *hist_sig_lep_rap = new TH1F("sig_lep_rap","sig_lep_rap", 40, -4.5 , 4.5);



  TH1F* h_cos_theta_star_lepminus = new TH1F("cos_theta_star_lepminus", "cos_theta_star_lepminus", 20, -1.0, 1.0); 
  TH1F* h_cos_theta_star_lepminus_cut1 = new TH1F("cos_theta_star_lepminus_cut1", "cos_theta_star_lepminus_cut1", 20, -1.0, 1.0); 


  TH2F* h2_cos_lep_PM = new TH2F("cos_lep_PM", "cos_lep_PM", 20,-1, 1, 20, -1, 1);



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
   TFile *outfile = new TFile("/data/SOM_Test/ALP_test/OutputFile_ALP_Z_125_Delphes.root","RECREATE");

  TTree *tree=new TTree("Tree","Signal");
  Float_t  Z_reco_pt = 0., LepP_pt=0., LepM_pt=0., MET_et = 0.0, Dphi_LepP_LepM=0., Deta_LepP_LepM=0., Dphi_Z_Diphoton = 0., Deta_Z_Diphoton=0.;
  Float_t LP=0., Frac_energy_diff=0., Cos_theta_star =0.0, Phi_star=0.;
  Int_t Total_Events=0 ;

  tree->Branch("Total_Events" ,     & Total_Events ,      "Total_Events/I");
  tree->Branch("Z_reco_pt",&Z_reco_pt,"Z_reco_pt/F");
  tree->Branch("LepP_pt",&LepP_pt,"LepP_pt/F");
  tree->Branch("LepM_pt",&LepM_pt,"LepM_pt/F");
  tree->Branch("MET_et",&MET_et,"MET_et/F");
  tree->Branch("Dphi_LepP_LepM",&Dphi_LepP_LepM,"Dphi_LepP_LepM/F");
  tree->Branch("Deta_LepP_LepM",&Deta_LepP_LepM,"Deta_LepP_LepM/F");
  tree->Branch("Dphi_Z_Diphoton",&Dphi_Z_Diphoton,"Dphi_Z_Diphoton/F");
  tree->Branch("Deta_Z_Diphoton",&Deta_Z_Diphoton,"Deta_Z_Diphoton/F");
  tree->Branch("LP",&LP,"LP/F");
  tree->Branch("Frac_energy_diff",&Frac_energy_diff,"Frac_energy_diff/F");
  tree->Branch("Cos_theta_star",&Cos_theta_star,"Cos_theta_star/F");
  tree->Branch("Phi_star",&Phi_star,"Phi_star/F");


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

int Z_idx = -1;
int lepton1_idx = -1;
int lepton2_idx = -1;
int ALP_idx = -1;
int light_dau1 = -1;
int light_dau2 = -1;
int ALP_gam1 = -1;
int ALP_gam2 = -1;
////////////////////////////////////////////////////////////////////
for(int i = 0; i < branchParticle->GetEntriesFast(); ++i) {
    particle = (GenParticle*) branchParticle->At(i);

    if (abs (particle->PID) == 24 && abs(particle->Status) == 62){
        Z_idx = i;
        //cout<<particle->Status<<"  "<<entry<<endl;
    }
    if ( (abs(particle->PID) == 25 ||  abs(particle->PID) == 9000005 ||  abs(particle->PID) == 36) && abs(particle->Status) == 62) {  // ALP
                ALP_idx= i;
    }
    if (abs(particle->PID) == 23 && abs(particle->Status) == 62) {  // Z boson (after FSR)
                Z_idx = i;
                // Check daughters
                int d1 = particle->D1;
                int d2 = particle->D2;
                
                if (d1 > 0 && d2 > 0) {

                    auto* daughter1= (GenParticle*) branchParticle->At(d1);
                    auto* daughter2= (GenParticle*) branchParticle->At(d2);
                    
                    // Check for leptonic decay (e, mu, tau)
                    if ( (abs(daughter1->PID) == 11 || abs(daughter1->PID) == 13 || abs(daughter1->PID) == 15) && daughter1->Charge == 1.0){
                        lepton1_idx = d1;
                   
                    }if ( (abs(daughter2->PID) == 11 || abs(daughter2->PID) == 13 || abs(daughter2->PID) == 15) && daughter2->Charge == -1.0) {
                        lepton2_idx = d2;
                    
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
n1++;
//cout<<Z_idx<<"  "<<ALP_idx<<"  "<<lepton1_idx<<"  "<<lepton2_idx<<"  "<<ALP_gam1<<"   "<<ALP_gam2<<endl;
////////////////////////Gen Level Checks/////////////////////
GenParticle *Z = (GenParticle*) branchParticle->At(Z_idx);
GenParticle *ALP = (GenParticle*) branchParticle->At(ALP_idx);
GenParticle *lepton1 = (GenParticle*) branchParticle->At(lepton1_idx); 
GenParticle *lepton2 = (GenParticle*) branchParticle->At(lepton2_idx); 
GenParticle *gamma1 = (GenParticle*) branchParticle->At(ALP_gam1); 
GenParticle *gamma2 = (GenParticle*) branchParticle->At(ALP_gam2);
if(lepton2->PT > lepton1->PT) {
    std::swap(lepton1, lepton2);
}

TLorentzVector lepP_vec_org_gen, lepM_vec_org_gen;
if (lepton1->PT > 20 && abs(lepton1->Eta) < 3.0 && gamma1->PT > 10 && gamma2->PT > 10 && abs(gamma1->Eta) < 2.5 && 
abs(gamma2->Eta) < 2.5  && lepton2->PT > 10 && abs(lepton2->Eta) < 3.0){ 

    TLorentzVector Z_vec_gen(Z->Px, Z->Py, Z->Pz, Z->E);
    TLorentzVector lepton2_vec_gen(lepton2->Px, lepton2->Py, lepton2->Pz, lepton2->E);
    TLorentzVector lepton1_vec_gen(lepton1->Px, lepton1->Py, lepton1->Pz, lepton1->E);
    TLorentzVector gam1_vec_gen(gamma1->Px, gamma1->Py, gamma1->Pz, gamma1->E);
    TLorentzVector gam2_vec_gen(gamma2->Px, gamma2->Py, gamma2->Pz, gamma2->E);
    TLorentzVector diphoton_gen = gam1_vec_gen + gam2_vec_gen ;
    if (lepton1_vec_gen.DeltaR(lepton2_vec_gen) < 0.3 || lepton1_vec_gen.DeltaR(gam1_vec_gen) < 0.3 || gam2_vec_gen.DeltaR(gam1_vec_gen) < 0.3  ) continue; 

n2++;

    lepP_vec_org_gen = lepton1_vec_gen;
    lepM_vec_org_gen = lepton2_vec_gen;


    // Reconstruct Z
    TLorentzVector Z_reco_gen = lepP_vec_org_gen + lepM_vec_org_gen;


    //cout<<Z_vec_gen.Vect().Mag()<<"  "<<Z_reco_gen.Vect().Mag()<<endl;
    hist_W_mass_gen->Fill( Z_vec_gen.M());
    if(abs(Z_reco_gen.M() - 91.1876) > 10) continue;
    hist_W_mass_gen1->Fill(Z_reco_gen.M());
    hist_W_rap_gen->Fill(Z_reco_gen.Rapidity());
    hist_W_pt3->Fill(Z_vec_gen.Pt());

    hist_MET_gen->Fill(Genmet->MET);
    hist_lep_pt_gen->Fill(lepton1->PT);
    hist_lep_rap_gen->Fill(lepP_vec_org_gen.Eta());
    hist_nu_pt->Fill(lepton2->PT);


    hist_dphi_lep_MET_gen->Fill(lepP_vec_org_gen.DeltaPhi(lepM_vec_org_gen));
    hist_deta_lep_MET_gen->Fill(lepP_vec_org_gen.Eta() - lepM_vec_org_gen.Eta());
    hist_dphi_Wrecopt_digam_gen->Fill(Z_reco_gen.DeltaPhi(diphoton_gen));
    hist_deta_Wrecopt_digam_gen->Fill(Z_reco_gen.Eta() - diphoton_gen.Eta());

    hist_mag_3momenta_Wgen->Fill(Z_vec_gen.Vect().Mag());
    hist_W_pt_gen->Fill(Z_vec_gen.Pt());
    
    TVector3 boostZ = Z_vec_gen.BoostVector();
    lepton1_vec_gen.Boost(-boostZ);
    TVector3 lep_dir_Zrest_gen = lepton1_vec_gen.Vect().Unit();   // lepton1 in Z rest frame
    TVector3 Z_dir_lab_gen     = Z_vec_gen.Vect().Unit();     // Z direction in lab
    double cos_theta_star_gen = lep_dir_Zrest_gen.Dot(Z_dir_lab_gen);
    h_cos_theta_star_gen->Fill(cos_theta_star_gen);
    double costhetastar_equivalant = abs(lepton1->E - lepton2->E)/Z_vec_gen.Vect().Mag();

    h_cos_theta_star_new_gen->Fill(costhetastar_equivalant);
    
    /////LP 
    double LP_gen= (lepP_vec_org_gen.Px()*Z_vec_gen.Px() + lepP_vec_org_gen.Py()*Z_vec_gen.Py())
                  / (Z_vec_gen.Px()*Z_vec_gen.Px() + Z_vec_gen.Py()*Z_vec_gen.Py());

    hist_LP_gen->Fill(LP_gen);

    if (Z_vec_gen.Pt() > 100) {
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

//=============================================== Leptons==============================
std::vector<fastjet::PseudoJet> Lepton_vec;
std::vector<int> Lepton_charge;
std::vector<int> Lepton_flavor; // 11 = e, 13 = mu

for(Int_t i = 0; i < branchElectron->GetEntriesFast(); ++i) {
    Electron *electron = (Electron*) branchElectron->At(i);
    if(electron->PT > 10.0 && fabs(electron->Eta) < 3.0 ) {
        TLorentzVector p4 = electron->P4();
        Lepton_vec.emplace_back(p4.Px(), p4.Py(), p4.Pz(), p4.E());

        Lepton_charge.push_back(electron->Charge);
        Lepton_flavor.push_back(11);
      }
    }
for(Int_t i = 0; i < branchMuon->GetEntriesFast(); ++i) {
    Muon *muon = (Muon*) branchMuon->At(i);
    if(muon->PT > 10.0 && fabs(muon->Eta) < 3.0) {
        TLorentzVector p4 = muon->P4();
        Lepton_vec.emplace_back(p4.Px(), p4.Py(), p4.Pz(), p4.E());

        Lepton_charge.push_back(muon->Charge);
        Lepton_flavor.push_back(13);
      }
    }
//cout<<Lepton_vec.size()<<endl;
hist_lep_size->Fill(Lepton_vec.size());
n3++;
////////////////////////////////////Delphes Level Z selection///////////////////////////
double cos_theta_star =-5, phi_star = -5, cos_theta_star_lepminus=-5;    
const double MZ = 91.1876;
const double Z_window = 10.0;
const double DR_cut = 0.3;  /////////////////////// 0.3 for comparison

bool foundZ_photons = false;
double best_diff = 1e9;

TLorentzVector Z_vec;
TLorentzVector lep_minus, lep_plus;
TLorentzVector gam1_best, gam2_best;

if (Lepton_vec.size() >= 2 && Photon_vec.size() >= 2){

    // loop over photon pairs
    for (size_t g1 = 0; g1 < Photon_vec.size(); ++g1) {
        for (size_t g2 = g1 + 1; g2 < Photon_vec.size(); ++g2) {

            TLorentzVector gam1_vec(Photon_vec[g1].px(), Photon_vec[g1].py(),
                                   Photon_vec[g1].pz(), Photon_vec[g1].e());

            TLorentzVector gam2_vec(Photon_vec[g2].px(), Photon_vec[g2].py(),
                                   Photon_vec[g2].pz(), Photon_vec[g2].e());

            // γ–γ separation
            if (gam1_vec.DeltaR(gam2_vec) < DR_cut) continue;

            // loop over lepton pairs
            for (size_t i = 0; i < Lepton_vec.size(); ++i) {
                for (size_t j = i + 1; j < Lepton_vec.size(); ++j) {

                    // SF
                    if (Lepton_flavor[i] != Lepton_flavor[j]) continue;

                    // OS
                    if (Lepton_charge[i] * Lepton_charge[j] != -1) continue;

                    TLorentzVector l1(Lepton_vec[i].px(), Lepton_vec[i].py(),
                                      Lepton_vec[i].pz(), Lepton_vec[i].e());

                    TLorentzVector l2(Lepton_vec[j].px(), Lepton_vec[j].py(),
                                      Lepton_vec[j].pz(), Lepton_vec[j].e());

                    // ℓ–γ separation
                    if (l1.DeltaR(gam1_vec) < DR_cut) continue;
                    if (l1.DeltaR(gam2_vec) < DR_cut) continue;
                    if (l2.DeltaR(gam1_vec) < DR_cut) continue;
                    if (l2.DeltaR(gam2_vec) < DR_cut) continue;

                    // pT cuts
                    double pt1 = l1.Pt();
                    double pt2 = l2.Pt();

                    if (std::max(pt1, pt2) < 20.0) continue;
                    if (std::min(pt1, pt2) < 10.0) continue;

                    // dilepton mass
                    TLorentzVector dilep = l1 + l2;
                    double mll = dilep.M();
                    double diff = fabs(mll - MZ);

                    // choose best Z candidate
                    if (diff < Z_window && diff < best_diff) {

                        best_diff = diff;
                        foundZ_photons = true;

                        Z_vec = dilep;
                        gam1_best = gam1_vec;
                        gam2_best = gam2_vec;

                        if (Lepton_charge[i] < 0){
                            lep_minus = l1;
                            lep_plus  = l2;
                        } else {
                            lep_minus = l2;
                            lep_plus  = l1;
                        }
                    }
                }
            }
        }
    }
}
if(foundZ_photons){ //////////////////////////Basic Selection
n4++;
    TLorentzVector diphoton = gam1_best + gam2_best;

    TLorentzVector Z_reco = lep_minus + lep_plus;

    //if (diphoton.DeltaR(Z_reco) < 0.4) continue; ////////////// cut on delta 


    hist_W_mass->Fill(Z_reco.M());
    hist_W_rap->Fill(Z_reco.Rapidity());

    hist_MET->Fill(met->MET);
    hist_DR_gen_Del_lep->Fill(lep_plus.DeltaR(lepP_vec_org_gen));


    hist_W_pt1->Fill(diphoton.Pt());

    // Reconstruct Z
    hist_dphi_lep_MET->Fill(abs(lep_minus.DeltaPhi(lep_plus)));

    hist_mag_3momenta_W->Fill(Z_reco.Vect().Mag());

    //cout<<Z->PT<<"   "<<Z_reco.Pt()<<"   "<<Z_reco.M()<<endl;

    Z_reco_pt = 0., LepP_pt = 0., MET_et=0., Dphi_LepP_LepM=0., Deta_LepP_LepM=0., Dphi_Z_Diphoton = 0., Deta_Z_Diphoton=0.;
    Z_reco_pt = Z_reco.Pt();
    LepP_pt = lep_plus.Pt();
    LepM_pt = lep_minus.Pt();
    MET_et = met->MET;
    Dphi_LepP_LepM= lep_plus.DeltaPhi(lep_minus);
    Deta_LepP_LepM = lep_plus.Eta() - lep_minus.Eta();

    Dphi_Z_Diphoton = Z_reco.DeltaPhi(diphoton);
    Deta_Z_Diphoton = Z_reco.Eta() - diphoton.Eta(); 

    hist_dphi_Wrecopt_digam->Fill(Dphi_Z_Diphoton);
    hist_deta_Wrecopt_digam->Fill(Deta_Z_Diphoton);

    //hist_dphi_Wrecopt_W->Fill(Z_reco.DeltaPhi(W_vec));
    h2_ptlep_vs_Wreco->Fill(lep_plus.Pt(), Z_reco.Pt());
    h2_ptlep_vs_digam->Fill(lep_plus.Pt(), diphoton.Pt());

    /////////Costheta*   
    //TLorentzVector W_vec(Z->Px, Z->Py, Z->Pz, Z->E);  
    //TVector3 boostW = W_vec.BoostVector();
    TLorentzVector lep_plus_CM = lep_plus;

    TVector3 boostW = Z_reco.BoostVector();
    lep_plus_CM.Boost(-boostW);

    TVector3 lep_dir_Zrest = lep_plus_CM.Vect().Unit();   // lepton1 in Z rest frame
    TVector3 Z_dir_lab     = Z_reco.Vect().Unit();     // Z direction in lab

    cos_theta_star = lep_dir_Zrest.Dot(Z_dir_lab);

    Cos_theta_star =0, Phi_star = 0;
    Cos_theta_star = cos_theta_star;
    h_cos_theta_star->Fill(cos_theta_star);    
    h2_cos_vs_ptW->Fill(Z_reco.Pt(), cos_theta_star);
    h2_cos_vs_ptlep->Fill(lep_plus.Pt(), cos_theta_star);


    // For lep minus costheta* 

    TLorentzVector lep_minus_CM = lep_minus;
    lep_minus_CM.Boost(-boostW);
    TVector3 lep_minus_dir_Zrest = lep_minus_CM.Vect().Unit();   // leptonminus in Z rest frame
    cos_theta_star_lepminus = lep_minus_dir_Zrest.Dot(Z_dir_lab);
    h_cos_theta_star_lepminus->Fill(cos_theta_star_lepminus);    
    //cout<<cos_theta_star<<"   "<<cos_theta_star_lepminus<<endl;
    h2_cos_lep_PM->Fill(cos_theta_star, cos_theta_star_lepminus);




////////////Costeta 1 
    // CM frame of final state objetcs
    TLorentzVector ALL_final = Z_reco + diphoton;
    TVector3 boostCM = ALL_final.BoostVector();
    // To get the Z dir in the CM frame
    TLorentzVector Z_reco_CM = Z_reco ;
    Z_reco_CM.Boost(-boostCM);
    // Dir of Z in the CM frame
    TVector3 Z_dir_CM = Z_reco_CM.Vect().Unit();
    double cos_theta1 = Z_dir_CM.CosTheta();  // wrt z-axis
    // Dir of Higgs in the CM frame
    TLorentzVector diphoton_CM = diphoton;
    diphoton_CM.Boost(-boostCM);
    TVector3 H_dir_CM = diphoton_CM.Vect().Unit();
    h_cos_theta1->Fill(cos_theta1);

    // check 
    //cout<<" In CM this should be -1 => "<<Z_dir_CM.Dot(H_dir_CM)<<endl;

//////Phi*
//phi_star azimuthal angle of the charged lepton1 in the Z rest frame
//measured with respect to the production plane
TVector3 beam_dir(0.0, 0.0, 1.0);   // proton beam direction
TVector3 y_hat = beam_dir.Cross(Z_dir_lab).Unit();
TVector3 x_hat = y_hat.Cross(Z_dir_lab).Unit();

phi_star = abs(atan2(lep_dir_Zrest.Dot(y_hat),lep_dir_Zrest.Dot(x_hat)));
Phi_star = phi_star;
h_phi_star_cut->Fill(Phi_star);

hist_lep_pt->Fill(lep_plus.Pt());
hist_sig_lep_rap->Fill(lep_plus.Eta());

hist_dphi_lep_gamma1->Fill(lep_plus.DeltaPhi(gam1_best));
hist_dphi_lep_gamma2->Fill(lep_plus.DeltaPhi(gam2_best));
hist_dR_lep_gamma1->Fill(lep_plus.DeltaR(gam1_best));

//////LP
LP=0.;
TVector3 pT_3lep = lep_plus.Vect();
TVector3 pT_3W   = Z_reco.Vect();

//LP= pT_3lep.Dot(pT_3W)/pT_3W.Mag2();
LP = (lep_plus.Px()*Z_reco.Px() + lep_plus.Py()*Z_reco.Py())
                  / (Z_reco.Px()*Z_reco.Px() + Z_reco.Py()*Z_reco.Py());

Frac_energy_diff =0;
Frac_energy_diff =abs(lep_plus.E() - lep_minus.E())/pT_3W.Mag() ;
h_cos_theta_star_new->Fill(Frac_energy_diff);

//cout<<cos_theta_star<<"   "<<2*(LP - 0.5)<<"  "<< abs(lep_plus.E() - MET_et)/Z_reco.Pt()<<endl;

hist_LP->Fill(LP);
h_cos_theta_star_from_LP->Fill(2*(LP - 0.5));
//if (Z_reco.Pt() > 50) cout<<Cos_theta_star<<"   "<<2*(LP - 0.5)<<"   "<<entry<<endl;

hist_ratio_ptl_ptW->Fill(lep_plus.Pt()/Z_reco.Pt());
hist_ratio_ptl_ptdigam->Fill(lep_plus.Pt()/diphoton.Pt());

    if(Z_reco.Pt() > 100){
    //cout<<Z_reco.M()<<endl;
    h_cos_theta_star_from_LP_cut1->Fill(2*(LP - 0.5));
    h_cos_theta_star_cut1->Fill(Cos_theta_star);
    hist_LP_cut1->Fill(LP);
    hist_W_pt3_cut1->Fill(Z_reco.Pt()); /////// var 
    hist_lep_pt_cut1->Fill(lep_plus.Pt());
    hist_ratio_ptl_ptW_cut1->Fill(lep_plus.Pt()/Z_reco.Pt());
    hist_dphi_lep_MET_cut1->Fill(abs(lep_plus.DeltaPhi(lep_minus)));
    h_cos_theta_star_new_cut1->Fill(abs(lep_plus.E() - lep_minus.E())/pT_3W.Mag());
    
    h_cos_theta_star_lepminus_cut1->Fill(cos_theta_star_lepminus);    
    
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
 TFile *fout = TFile::Open("/data/SOM_Test/ALP_test/Result_ALP_Z_125_Delphes.root","RECREATE");

fout->cd();
hist_nu_pz_gen->Write();
hist_W_rap->Write();
hist_W_rap_gen->Write();
hist_W_rap_gen1->Write();
h2_cos_vs_ptW->Write();
h2_cos_vs_ptlep->Write();

hist_W_mass_gen->Write();
hist_W_mass_gen1->Write();
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
hist_deta_Wrecopt_digam_gen->Write();
hist_dphi_lep_MET_gen->Write();
hist_deta_lep_MET_gen->Write();
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


h2_cos_lep_PM->Write();
h_cos_theta_star_lepminus->Write();
h_cos_theta_star_lepminus_cut1->Write();


hist_cut_flow->Write();
delete fout;
}// End of AnalyseEvents

void ALP_AP_Z(const char *inputFile)
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