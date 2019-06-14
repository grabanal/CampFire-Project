//
//  ChargedHiggsMC_code.cpp
//  
//
//  Created by Utsav Patel on 4/22/19.
//

#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include "TLatex.h"
#include <algorithm>
#include "math.h"
#include "TRandom.h"
#include "time.h"
using namespace std;

int maximum(float,float,float);
int minimum(float,float,float);
int jet_cuts(float,float,float,float,float,float);
//int * photon_filter(int,int,int);

void ChargedHiggsMC_code(){
    
    double pi = 3.141592653589793;
    
    clock_t code_start_time = clock();
    
    std::ifstream ifs("root_files.txt");
    std::string line;
    
    while(std::getline(ifs, line)){
        clock_t start_time = clock();
        printf("This is the analysed file: %s\n",line.c_str());
        
        TFile *myROOTFile = TFile::Open(line.c_str());
        
        TH1F *PtLGammaNu = new TH1F("PtLGammaNu", "P_{T}^{l+#gamma+MET}", 60, 0, 600);
        PtLGammaNu->SetTitle("P_{T}^{l+#gamma+MET};GeV/c;Event Fraction/10GeV");
        
        TH1F *MET = new TH1F("MET", "Missing E_{T}", 40, 0, 400);
        MET->SetTitle("Missing E_{T};GeV;Event Fraction/10GeV");
        
        TH1F *HT = new TH1F("HT", "H_{T}", 70, 0, 700);
        HT->SetTitle("H_{T};GeV/c;Event Fraction/10GeV");
        
        TH1F *Plq = new TH1F("Plq", "P_{l}.q", 100, 0, 10000);
        Plq->SetTitle("P_{l} #dot q;(GeV/c)^{2};Event Fraction/10GeV");

//        These are general histograms that fit every signal
        TH1F *ParticleCounter = new TH1F("ParticleCounter", "Number of Particles", 5, 0, 5);
        ParticleCounter->SetTitle("Number of Particles;Particle Type;Number of Events");
        TH1F *LeptonPt = new TH1F("LeptonPt", "P_{T} of the Lepton", 40, 0, 400);
        LeptonPt->SetTitle("P_{T} of the Lepton;GeV/c;Event Fraction/10GeV/c");
        TH1F *LeptonP = new TH1F("LeptonP", "Momentum Distribution of the Lepton", 60, 0, 600);
        LeptonP->SetTitle("Momentum Distribution of the Lepton;GeV/c;Event Fraction/10GeV/c");
        TH1F *PhotonEt = new TH1F("PhotonEt", "E_{T} of the Photon", 40, 0, 400);
        PhotonEt->SetTitle("E_{T} of the Photon;GeV;Event Fraction/10GeV");
        TH1F *PhotonE = new TH1F("PhotonE", "E Distribution of the Photon", 60, 0, 600);
        PhotonE->SetTitle("E Distribution of the Photon;GeV;Event Fraction/10GeV");
        TH1F *MtLPN = new TH1F("MtLPN", "Transverse Mass of the Lepton, Photon and Neutrino", 40, 0, 400);
        MtLPN->SetTitle("Transverse Mass of the Lepton, Photon and Neutrino;GeV/c^2;Event Fraction/10GeV/c^2");
        TH1F *MtLN = new TH1F("MtLN", "Transverse Mass of the Lepton and Neutrino", 50, 0, 200);
        MtLN->SetTitle("Transverse Mass of the Lepton and Neutrino;GeV/c^2;Event Fraction/4GeV/c^2");
        TH1F *InvMassLPN = new TH1F("InvMassLPN", "Invariant Mass of the Lepton, Photon and Neutrino", 50, 0, 200);
        InvMassLPN->SetTitle("Invariant Mass of the Lepton, Photon and Neutrino;GeV/c^2;Event Fraction/4GeV/c^2");
        TH1F *InvMassLN = new TH1F("InvMassLN", "Invariant Mass of the Lepton and Neutrino", 60, 0, 120);
        InvMassLN->SetTitle("Invariant Mass of the Lepton and Neutrino;GeV/c^2;Event Fraction/2GeV/c^2");
        
        TH1F *MissingE = new TH1F("MissingE", "Missing Energy", 40, 0, 400);
        MissingE->SetTitle("Missing Energy;GeV;Event Fraction/10GeV");
        TH1F *MEwLepton = new TH1F("MEwLepton", "Missing Energy Associated with Leptonic Decay", 60, 0, 600);
        MEwLepton->SetTitle("Missing Energy Associated with Leptonic Decay;GeV;Event Fraction/10GeV");
        
        TH1F *DeltaPhiLP = new TH1F("DeltaPhiLP", "#Delta#phi Between Lepton and Photon", 50, 0, 3.2);
        DeltaPhiLP->SetTitle("#Delta#phi Between Lepton and Photon;Radians;Event Fraction");
        TH1F *DeltaPhiLN = new TH1F("DeltaPhiLN", "#Delta#phi Between Lepton and Neutrino", 50, 0, 3.2);
        DeltaPhiLN->SetTitle("#Delta#phi Between Lepton and Neutrino;Radians;Event Fraction");
        TH1F *DeltaEtaLP = new TH1F("DeltaEtaLP", "#Delta#eta Between Lepton and Photon", 50, 0, 5);
        DeltaEtaLP->SetTitle("#Delta#eta Between Lepton and Photon;;Event Fraction");
        TH1F *DeltaEtaLN = new TH1F("DeltaEtaLN", "#Delta#eta Between Lepton and Neutrino", 50, 0, 5);
        DeltaEtaLN->SetTitle("#Delta#eta Between Lepton and Neutrino;;Event Fraction");
        TH1F *DeltaRLP = new TH1F("DeltaRLP", "#Delta R Between Lepton and Photon", 40, 0, 4);
        DeltaRLP->SetTitle("#Delta R Between Lepton and Photon;;Event Fraction");
        TH1F *DeltaRLN = new TH1F("DeltaRLN", "#Delta R Between Lepton and Neutrino", 40, 0, 4);
        DeltaRLN->SetTitle("#Delta R Between Lepton and Neutrino;;Event Fraction");
        
        TH1F *LeptonPhi = new TH1F("LeptonPhi", "#phi Distribution of the Lepton", 50, -3.2, 3.2);
        LeptonPhi->SetTitle("#phi Distribution of the Lepton;Radians;Event Fraction");
        TH1F *LeptonEta = new TH1F("LeptonEta", "#eta Distribution of the Lepton", 50, 0, 2.5);
        LeptonEta->SetTitle("#eta Distribution of the Lepton;;Event Fraction");
        TH1F *PhotonPhi = new TH1F("PhotonPhi", "#phi Distribution of the Emitted Photon(s)", 50, -3.2, 3.2);
        PhotonPhi->SetTitle("#phi Distribution of the Emitted Photon(s);Radians;Event Fraction");
        TH1F *PhotonEta = new TH1F("PhotonEta", "#eta Distribution of the Emitted Photon(s)", 50, 0, 2.5);
        PhotonEta->SetTitle("#eta Distribution of the Emitted Photon(s);;Event Fraction");
        
//        These are histograms that only apply to the first signal (charged and uncharged higgs)
        TH1F *DeltaPhi2Photons = new TH1F("DeltaPhi2Photons", "#Delta#phi Between Photons from H^{0}", 50, 0, 3.2);
        DeltaPhi2Photons->SetTitle("#Delta#phi Between Photons from H^{0};;Event Fraction");
        TH1F *DeltaEta2Photons = new TH1F("DeltaEta2Photons", "#Delta#eta Between Photons from H^{0}", 50, 0, 5);
        DeltaEta2Photons->SetTitle("#Delta#eta Between Photons from H^{0};;Event Fraction");
        TH1F *DeltaR2Photons = new TH1F("DeltaR2Photons", "#Delta R Between Photons from H^{0}", 40, 0, 4);
        DeltaR2Photons->SetTitle("#Delta R Between Photons from H^{0};;Event Fraction");
        TH1F *ClosestPhoton = new TH1F("ClosestPhoton", "Closest Photon to the Lepton", 3, 0, 3);
        ClosestPhoton->SetTitle("Closest Photon to the Lepton;;Number of Events");
        TH1F *InvMass2Photons = new TH1F("InvMass2Photons", "Invariant Mass of the Photons from H^{0}", 50, 0, 200);
        InvMass2Photons->SetTitle("Invariant Mass of the Photons from H^{0};GeV/c^2;Event Fraction/4GeV/c^2");
        TH1F *Mt2Photons = new TH1F("Mt2Photons", "Transverse Mass of the Photons from H^{0}", 50, 0, 200);
        Mt2Photons->SetTitle("Transverse Mass of the Photons from H^{0};GeV/c^2;Event Fraction/4GeV/c^2");
        
//        These are histograms that only apply to the second signal (charged and doubly charged higgs)
        TH1F *DeltaPhiJets1 = new TH1F("DeltaPhiJets1", "#Delta#phi Between Jets from the First W Decay", 50, 0, 3.2);
        DeltaPhiJets1->SetTitle("#Delta#phi Between Jets from the First W Decay;;Event Fraction");
        TH1F *DeltaPhiJets2 = new TH1F("DeltaPhiJets2", "#Delta#phi Between Jets from the Second W Decay", 50, 0, 3.2);
        DeltaPhiJets2->SetTitle("#Delta#phi Between Jets from the Second W Decay;;Event Fraction");
        TH1F *DeltaEtaJets1 = new TH1F("DeltaEtaJets1", "#Delta#eta Between Jets from the First W Decay", 50, 0, 5);
        DeltaEtaJets1->SetTitle("#Delta#eta Between Jets from the First W Decay;;Event Fraction");
        TH1F *DeltaEtaJets2 = new TH1F("DeltaEtaJets2", "#Delta#eta Between Jets from the Second W Decay", 50, 0, 5);
        DeltaEtaJets2->SetTitle("#Delta#eta Between Jets from the Second W Decay;;Event Fraction");
        TH1F *DeltaRJets1 = new TH1F("DeltaRJets1", "#Delta R Between Jets from the First W Decay", 40, 0, 4);
        DeltaRJets1->SetTitle("#Delta R Between Jets from the First W Decay;;Event Fraction");
        TH1F *DeltaRJets2 = new TH1F("DeltaRJets2", "#Delta R Between Jets from the Second W Decay", 40, 0, 4);
        DeltaRJets2->SetTitle("#Delta R Between Jets from the Second W Decay;;Event Fraction");
        TH1F *InvMassJets1 = new TH1F("InvMassJets1", "Invariant Mass of the Jets from the First W Decay", 50, 0, 100);
        InvMassJets1->SetTitle("Invariant Mass of the Jets from the First W Decay;Event Fraction/2GeV/c^2");
        TH1F *InvMassJets2 = new TH1F("InvMassJets2", "Invariant Mass of the Jets from the Second W Decay", 50, 0, 100);
        InvMassJets2->SetTitle("Invariant Mass of the Jets from the Second W Decay;Event Fraction/2GeV/c^2");
        TH1F *MtJets1 = new TH1F("MtJets1", "Transverse Mass of the Jets from the First W Decay", 50, 0, 100);
        MtJets1->SetTitle("Transverse Mass of the Jets from the First W Decay;GeV/c^2;Event Fraction/2GeV/c^2");
        TH1F *MtJets2 = new TH1F("MtJets2", "Transverse Mass of the Jets from the Second W Decay", 50, 0, 100);
        MtJets2->SetTitle("Transverse Mass of the Jets from the Second W Decay;GeV/c^2;Event Fraction/2GeV/c^2");
        
//        These are histograms that only apply to the third signal (two charged higgs)
        TH1F *DeltaPhiSepPhotons = new TH1F("DeltaPhiSepPhotons", "#Delta#phi Between Photons from Both H^{#pm} Decays", 50, 0, 3.2);
        DeltaPhiSepPhotons->SetTitle("#Delta#phi Between Photons from Both H^{#pm} Decays;;Event Fraction");
        TH1F *DeltaEtaSepPhotons = new TH1F("DeltaEtaSepPhotons", "#Delta#eta Between Photons from Both H^{#pm} Decays", 50, 0, 5);
        DeltaEtaSepPhotons->SetTitle("#Delta#eta Between Photons from Both H^{#pm} Decays;;Event Fraction");
        TH1F *DeltaRSepPhotons = new TH1F("DeltaRSepPhotons", "#Delta R Between Photons from Both H^{#pm} Decays", 40, 0, 4);
        DeltaRSepPhotons->SetTitle("#Delta R Between Photons from Both H^{#pm} Decays;;Event Fraction");
        
//        Extra histograms
        TH1F *Lepton1Pt = new TH1F("Lepton1Pt", "P_{T} of Lepton from the Higgs", 40, 0, 400);
        Lepton1Pt->SetTitle("P_{T} of Lepton from the Higgs;GeV/c;Event Fraction/10GeV/c");
        TH1F *Lepton2Pt = new TH1F("Lepton2Pt", "P_{T} of Lepton from the Higgs", 40, 0, 400);
        Lepton2Pt->SetTitle("P_{T} of Lepton from the Higgs;GeV/c;Event Fraction/10GeV/c");
        TH1F *Photon1E = new TH1F("Photon1E", "E of Photon from the Higgs", 40, 0, 400);
        Photon1E->SetTitle("E of Photon from the Higgs;GeV/c;Event Fraction/10GeV/c");
        TH1F *Photon2E = new TH1F("Photon2E", "E of Photon from the Higgs", 40, 0, 400);
        Photon2E->SetTitle("E of Photon from the Higgs;GeV/c;Event Fraction/10GeV/c");
        
        TH1F *HighestPtLepton = new TH1F("HighestPtLepton", "Lepton with the Highest P_{T}", 3, 0, 3);
        HighestPtLepton->SetTitle("Lepton with the Highest P_{T};Lepton Number;Fraction of Total");
        TH1F *ClosestNeutrino = new TH1F("ClosestNeutrino", "Neutrino Closest to the Highest P_{T} Lepton", 3, 0, 3);
        ClosestNeutrino->SetTitle("Neutrino Closest to the Highest P_{T} Lepton;Neutrino Number;Fraction of Total");
        TH1F *SameNeutrino = new TH1F("SameNeutrino", "Same Neutrino", 2, 0, 2);
        SameNeutrino->SetTitle("Same Neutrino;Yes/No;Event Fraction");
        
        
        TTree *particle_Tree = (TTree*)myROOTFile->Get("LHEF");
        
        TLeaf *pID = particle_Tree->GetLeaf("Particle/Particle.PID");
        TLeaf *pt = particle_Tree->GetLeaf("Particle/Particle.PT");
        TLeaf *mass = particle_Tree->GetLeaf("Particle/Particle.M");
        TLeaf *phi = particle_Tree->GetLeaf("Particle/Particle.Phi");
        TLeaf *energy = particle_Tree->GetLeaf("Particle/Particle.E");
        TLeaf *px = particle_Tree->GetLeaf("Particle/Particle.Px");
        TLeaf *py = particle_Tree->GetLeaf("Particle/Particle.Py");
        TLeaf *pz = particle_Tree->GetLeaf("Particle/Particle.Pz");
        TLeaf *eta = particle_Tree->GetLeaf("Particle/Particle.Eta");
        TLeaf *mother1 = particle_Tree->GetLeaf("Particle/Particle.Mother1");
        TLeaf *mother2 = particle_Tree->GetLeaf("Particle/Particle.Mother2");
        TRandom *rand_num = new TRandom();
        int num_events = particle_Tree->GetEntries();
        int num_cuts[5] = {0};
        int num_fid_cuts[5] = {0};
        
        float hcal_eta = 4.9;
        float tracker_eta = 2.5;
        float ecal_eta = 3.2;
        
        int num_all_cuts = 0;
        for(int i=0;i<num_events;i++){
            cout << i << '\r' << flush;
            //Array For Holding Information
            //Note for these arrays, the order is px, py, pz, pt, eta, phi, energy, mass
            float lepton1[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Lepton 1
            float lepton2[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Lepton 2
            float lepton3[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Lepton 3
            float neutrino1[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Neutrino 1
            float neutrino2[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Neutrino 2
            float neutrino3[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Neutrino 3
            float quark1[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Quark (not bottom) if present.
            float quark2[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Second Quark (will accompany first quark)
            float quark3[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Third Quark if present
            float quark4[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Fourth Quark (will accompany third quark)
            float photon1[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //First Photon
            float photon2[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Second Photon if present
            float photon3[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Third Photon if present
            float bottom1[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Bottom quark
            float bottom2[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //Anti-bottom quark
            
            float lepton[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //This is the lepton with the highest pt
            float neutrino[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //This is the associated neutrino
            float photon[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //This is the closest photon to the lepton with the highest pt
            
            bool lepton_present[2] = {false,false};
            bool neutrino_present[3] = {false,false,false};
            bool quark_present[3] = {false,false,false};
            bool bottom_present = false;
            bool photon_present[3] = {false,false,false};
            bool photon_passed[3] = {false,false,false};

            bool made_cuts[4] = {false,false,false,false};
            
            particle_Tree->GetEntry(i);
            int num_particles = pID->GetLen();
            for(int j=0;j<num_particles;j++){
                int particleID = pID->GetValue(j);
                if(particleID==11 || particleID==-11 || particleID==13 || particleID==-13 || particleID==15 || particleID==-15){
                    if(not lepton_present[0]){
                        lepton_present[0] = true;
                        lepton1[0] = px->GetValue(j);
                        lepton1[1] = py->GetValue(j);
                        lepton1[2] = pz->GetValue(j);
                        lepton1[3] = pt->GetValue(j);
                        lepton1[4] = energy->GetValue(j);
                        lepton1[5] = mass->GetValue(j);
                        lepton1[6] = phi->GetValue(j);
                        lepton1[7] = eta->GetValue(j);
                    }else if(not lepton_present[1]){
                        lepton_present[1] = true;
                        lepton2[0] = px->GetValue(j);
                        lepton2[1] = py->GetValue(j);
                        lepton2[2] = pz->GetValue(j);
                        lepton2[3] = pt->GetValue(j);
                        lepton2[4] = energy->GetValue(j);
                        lepton2[5] = mass->GetValue(j);
                        lepton2[6] = phi->GetValue(j);
                        lepton2[7] = eta->GetValue(j);
                    }else{
                        lepton3[0] = px->GetValue(j);
                        lepton3[1] = py->GetValue(j);
                        lepton3[2] = pz->GetValue(j);
                        lepton3[3] = pt->GetValue(j);
                        lepton3[4] = energy->GetValue(j);
                        lepton3[5] = mass->GetValue(j);
                        lepton3[6] = phi->GetValue(j);
                        lepton3[7] = eta->GetValue(j);
                    }
                    if(particleID==11 || particleID==-11){
                        ParticleCounter->Fill(0);
                    }else if(particleID==13 || particleID==-13){
                        ParticleCounter->Fill(1);
                    }else{
                        ParticleCounter->Fill(2);
                    }
                }else if(particleID==12 || particleID==-12 || particleID==14 || particleID==-14 || particleID==16 || particleID==-16){
                    if(not neutrino_present[0]){
                        neutrino_present[0] = true;
                        neutrino1[0] = px->GetValue(j);
                        neutrino1[1] = py->GetValue(j);
                        neutrino1[2] = pz->GetValue(j);
                        neutrino1[3] = pt->GetValue(j);
                        neutrino1[4] = energy->GetValue(j);
                        neutrino1[5] = mass->GetValue(j);
                        neutrino1[6] = phi->GetValue(j);
                        neutrino1[7] = eta->GetValue(j);
                    }else if(not neutrino_present[1]){
                        neutrino_present[1] = true;
                        neutrino2[0] = px->GetValue(j);
                        neutrino2[1] = py->GetValue(j);
                        neutrino2[2] = pz->GetValue(j);
                        neutrino2[3] = pt->GetValue(j);
                        neutrino2[4] = energy->GetValue(j);
                        neutrino2[5] = mass->GetValue(j);
                        neutrino2[6] = phi->GetValue(j);
                        neutrino2[7] = eta->GetValue(j);
                    }else{
                        neutrino_present[2] = true;
                        neutrino3[0] = px->GetValue(j);
                        neutrino3[1] = py->GetValue(j);
                        neutrino3[2] = pz->GetValue(j);
                        neutrino3[3] = pt->GetValue(j);
                        neutrino3[4] = energy->GetValue(j);
                        neutrino3[5] = mass->GetValue(j);
                        neutrino3[6] = phi->GetValue(j);
                        neutrino3[7] = eta->GetValue(j);
                    }
                }else if(particleID==22){
                    ParticleCounter->Fill(3);
                    if(not photon_present[0]){
                        photon_present[0] = true;
                        photon1[0] = px->GetValue(j);
                        photon1[1] = py->GetValue(j);
                        photon1[2] = pz->GetValue(j);
                        photon1[3] = pt->GetValue(j);
                        photon1[4] = energy->GetValue(j);
                        photon1[5] = mass->GetValue(j);
                        photon1[6] = phi->GetValue(j);
                        photon1[7] = eta->GetValue(j);
                        if(photon1[3]>25.0 and fabs(photon1[7])<2.5){
                            photon_passed[0] = true;
                        }
                    }else if(not photon_present[1]){
                        photon_present[1] = true;
                        photon2[0] = px->GetValue(j);
                        photon2[1] = py->GetValue(j);
                        photon2[2] = pz->GetValue(j);
                        photon2[3] = pt->GetValue(j);
                        photon2[4] = energy->GetValue(j);
                        photon2[5] = mass->GetValue(j);
                        photon2[6] = phi->GetValue(j);
                        photon2[7] = eta->GetValue(j);
                        if(photon2[3]>25.0 and fabs(photon2[7])<2.5){
                            photon_passed[1] = true;
                        }
                    }else{
                        photon_present[2] = true;
                        photon3[0] = px->GetValue(j);
                        photon3[1] = py->GetValue(j);
                        photon3[2] = pz->GetValue(j);
                        photon3[3] = pt->GetValue(j);
                        photon3[4] = energy->GetValue(j);
                        photon3[5] = mass->GetValue(j);
                        photon3[6] = phi->GetValue(j);
                        photon3[7] = eta->GetValue(j);
                        if(photon3[3]>25.0 and fabs(photon3[7])<2.5){
                            photon_passed[2] = true;
                        }
                    }
                }else if(particleID==5 || particleID==-5){
                    ParticleCounter->Fill(4);
                    if(not bottom_present){
                        bottom_present = true;
                        bottom1[0] = px->GetValue(j);
                        bottom1[1] = py->GetValue(j);
                        bottom1[2] = pz->GetValue(j);
                        bottom1[3] = pt->GetValue(j);
                        bottom1[4] = energy->GetValue(j);
                        bottom1[5] = mass->GetValue(j);
                        bottom1[6] = phi->GetValue(j);
                        bottom1[7] = eta->GetValue(j);
                    }else{
                        bottom2[0] = px->GetValue(j);
                        bottom2[1] = py->GetValue(j);
                        bottom2[2] = pz->GetValue(j);
                        bottom2[3] = pt->GetValue(j);
                        bottom2[4] = energy->GetValue(j);
                        bottom2[5] = mass->GetValue(j);
                        bottom2[6] = phi->GetValue(j);
                        bottom2[7] = eta->GetValue(j);
                    }
                }else if((particleID==1 || particleID==2 || particleID==3 || particleID==4 || particleID==-1 || particleID==-2 || particleID==-3 || particleID==-4) and j>1){
                    ParticleCounter->Fill(4);
                    if(not quark_present[0]){
                        quark_present[0] = true;
                        quark1[0] = px->GetValue(j);
                        quark1[1] = py->GetValue(j);
                        quark1[2] = pz->GetValue(j);
                        quark1[3] = pt->GetValue(j);
                        quark1[4] = energy->GetValue(j);
                        quark1[5] = mass->GetValue(j);
                        quark1[6] = phi->GetValue(j);
                        quark1[7] = eta->GetValue(j);
                    }else if(not quark_present[1]){
                        quark_present[1] = true;
                        quark2[0] = px->GetValue(j);
                        quark2[1] = py->GetValue(j);
                        quark2[2] = pz->GetValue(j);
                        quark2[3] = pt->GetValue(j);
                        quark2[4] = energy->GetValue(j);
                        quark2[5] = mass->GetValue(j);
                        quark2[6] = phi->GetValue(j);
                        quark2[7] = eta->GetValue(j);
                    }else if(not quark_present[2]){
                        quark_present[2] = true;
                        quark3[0] = px->GetValue(j);
                        quark3[1] = py->GetValue(j);
                        quark3[2] = pz->GetValue(j);
                        quark3[3] = pt->GetValue(j);
                        quark3[4] = energy->GetValue(j);
                        quark3[5] = mass->GetValue(j);
                        quark3[6] = phi->GetValue(j);
                        quark3[7] = eta->GetValue(j);
                    }else{
                        quark4[0] = px->GetValue(j);
                        quark4[1] = py->GetValue(j);
                        quark4[2] = pz->GetValue(j);
                        quark4[3] = pt->GetValue(j);
                        quark4[4] = energy->GetValue(j);
                        quark4[5] = mass->GetValue(j);
                        quark4[6] = phi->GetValue(j);
                        quark4[7] = eta->GetValue(j);
                    }
                }
            }
            
//            Selecting the lepton
            int lepton_highest_pt = maximum(lepton1[3],lepton2[3],lepton3[3]);
            if(lepton_highest_pt==1){
                for(int k=0;k<8;k++){
                    lepton[k]=lepton1[k];
                }
            }else if(lepton_highest_pt==2){
                for(int k=0;k<8;k++){
                    lepton[k]=lepton2[k];
                }
            }else{
                for(int k=0;k<8;k++){
                    lepton[k]=lepton3[k];
                }
            }
//            Selecting the photon
//            Have the photon pass and then select the closest one
            double deltar_photon[3] = {100.0,100.0,100.0};
            if(photon_present[0] and photon_passed[0]){
                double DeltaPhiLP1_val = fabs(lepton[6]-photon1[6]);
                if(DeltaPhiLP1_val > pi){
                    DeltaPhiLP1_val = 2*pi - DeltaPhiLP1_val;
                }
                deltar_photon[0] = sqrt(pow(DeltaPhiLP1_val,2)+pow((photon1[7]-lepton[7]),2));
            }
            if(photon_present[1] and photon_passed[1]){
                double DeltaPhiLP2_val = fabs(lepton[6]-photon2[6]);
                if(DeltaPhiLP2_val > pi){
                    DeltaPhiLP2_val = 2*pi - DeltaPhiLP2_val;
                }
                deltar_photon[1] = sqrt(pow(DeltaPhiLP2_val,2)+pow((photon2[7]-lepton[7]),2));
            }
            if(photon_present[2] and photon_passed[2]){
                double DeltaPhiLP3_val = fabs(lepton[6]-photon3[6]);
                if(DeltaPhiLP3_val > pi){
                    DeltaPhiLP3_val = 2*pi - DeltaPhiLP3_val;
                }
                deltar_photon[2] = sqrt(pow(DeltaPhiLP3_val,2)+pow((photon3[7]-lepton[7]),2));
            }
            int photon_lowest_deltar = minimum(deltar_photon[0],deltar_photon[1],deltar_photon[2]);
            
            if(photon_lowest_deltar==1 and photon_passed[0]){
                for(int k=0;k<8;k++){
                    photon[k]=photon1[k];
                }
            }else if(photon_lowest_deltar==2 and photon_passed[1]){
                for(int k=0;k<8;k++){
                    photon[k]=photon2[k];
                }
            }else if(photon_lowest_deltar==3 and photon_passed[2]){
                for(int k=0;k<8;k++){
                    photon[k]=photon3[k];
                }
            }
            
            
            
//            Select the closest neutrino to the lepton
            float deltar_neutrino[3] = {100.0,100.0,100.0};
            if(neutrino_present[0]){
                double DeltaPhiLN1_val = fabs(lepton[6]-neutrino1[6]);
                if(DeltaPhiLN1_val > pi){
                    DeltaPhiLN1_val = 2*pi - DeltaPhiLN1_val;
                }
                deltar_neutrino[0] = sqrt(pow(DeltaPhiLN1_val,2)+pow((neutrino1[7]-lepton[7]),2));
            }
            if(neutrino_present[1]){
                double DeltaPhiLN2_val = fabs(lepton[6]-neutrino2[6]);
                if(DeltaPhiLN2_val > pi){
                    DeltaPhiLN2_val = 2*pi - DeltaPhiLN2_val;
                }
                deltar_neutrino[1] = sqrt(pow(DeltaPhiLN2_val,2)+pow((neutrino2[7]-lepton[7]),2));
            }
            if(neutrino_present[2]){
                double DeltaPhiLN3_val = lepton[6]-neutrino3[6];
                if(DeltaPhiLN3_val > pi){
                    DeltaPhiLN3_val = 2*pi - DeltaPhiLN3_val;
                }
                deltar_neutrino[2] = sqrt(pow(DeltaPhiLN3_val,2)+pow((neutrino3[7]-lepton[7]),2));
            }
            
            int neutrino_lowest_deltar = minimum(deltar_neutrino[0],deltar_neutrino[1],deltar_neutrino[2]);
            
            
            if(neutrino_lowest_deltar==1){
                for(int k=0;k<8;k++){
                    neutrino[k]=neutrino1[k];
                }
            }else if(neutrino_lowest_deltar==2){
                for(int k=0;k<8;k++){
                    neutrino[k]=neutrino2[k];
                }
            }else if(neutrino_lowest_deltar==3){
                for(int k=0;k<8;k++){
                    neutrino[k]=neutrino3[k];
                }
            }
            
            num_cuts[0] += 1;
            
            float missing_et = sqrt(pow((neutrino1[0]+neutrino2[0]+neutrino3[0]),2) + pow((neutrino1[1]+neutrino2[1]+neutrino3[1]),2)); //Assume neutrino mass is approximately zero
            float h_t = lepton1[3]+lepton2[3]+lepton3[3]+bottom1[3]+bottom2[3]+quark1[3]+quark2[3]+quark3[3]+quark4[3]+photon1[3]+photon2[3]+photon3[3];
            float ptlgnu = sqrt(pow((lepton[0]+photon[0]),2) + pow((lepton[1]+photon[1]),2) + pow((neutrino1[0]+neutrino2[0]+neutrino3[0]),2) + pow((neutrino1[1]+neutrino2[1]+neutrino3[1]),2));
            float pl_q = lepton[4]*photon[4] - (lepton[0]*photon[0]+lepton[1]*photon[1]+lepton[2]*photon[2]);
            
//            First cut
            if(lepton[3]>25.0 and fabs(lepton[7])<2.5){
                num_cuts[1] += 1;
                made_cuts[0] = true;
            }
//            Second cut
            if(photon_passed[0] or photon_passed[1] or photon_passed[2]){
                num_cuts[2] += 1;
                made_cuts[1] = true;
            }
//            Third cut
            int num_qualifying_jets = jet_cuts(quark1[3],quark2[3],quark3[3],quark4[3],bottom1[3],bottom2[3]);
            if(num_qualifying_jets<=2){
                num_cuts[3] += 1;
                made_cuts[2] = true;
            }
//            Fourth cut
            if(not bottom_present){
                num_cuts[4] += 1;
                made_cuts[3] = true;
            }else if((rand_num->Rndm(i)>0.7) and (rand_num->Rndm(num_events*rand_num->Rndm(i))>0.7)){
                num_cuts[4] += 1;
                made_cuts[3] = true;
            }
            
            if((made_cuts[0] && made_cuts[1]) && (made_cuts[2] && made_cuts[3])){
                num_all_cuts += 1;
//                These are where I add every other plot. These can be moved to the fiducial cuts block if need be
                
            }
            
            
            if(fabs(photon[7]) < ecal_eta and fabs(lepton[7]) < tracker_eta and fabs(quark1[7]) < hcal_eta and fabs(quark2[7]) < hcal_eta and fabs(quark3[7]) < hcal_eta and fabs(quark4[7]) < hcal_eta and fabs(bottom1[7]) < hcal_eta and fabs(bottom2[7]) < hcal_eta){
                
                num_fid_cuts[0] += 1;
                
                if(made_cuts[0]){
                    num_fid_cuts[1] += 1;
                    if(made_cuts[1]){
                        num_fid_cuts[2] += 1;
                        if(made_cuts[2]){
                            num_fid_cuts[3] += 1;
                            if(made_cuts[3]){
                                num_fid_cuts[4] += 1;
                                PtLGammaNu->Fill(ptlgnu);
                                MET->Fill(missing_et);
                                HT->Fill(h_t);
                                Plq->Fill(pl_q);
                                LeptonPt->Fill(lepton[3]);
                                LeptonP->Fill(sqrt(pow(lepton[3],2)+pow(lepton[2],2)));
                                PhotonEt->Fill(photon[3]); //Photon Pt is equal to Photon Et since the photon is massless
                                PhotonE->Fill(photon[4]);
                                double et_lepton = sqrt(pow(lepton[5],2)+pow(lepton[0],2)+pow(lepton[1],2));
                                double et_photon = sqrt(pow(photon[5],2)+pow(photon[0],2)+pow(photon[1],2));
                                double et_neutrino = sqrt(pow(neutrino[5],2)+pow(neutrino[0],2)+pow(neutrino[1],2));
                                MtLPN->Fill(sqrt(pow((et_lepton+et_photon+et_neutrino),2) - pow((lepton[0]+photon[0]+neutrino[0]),2) - pow((lepton[1]+photon[1]+neutrino[1]),2)));
                                MtLN->Fill(sqrt(pow((et_lepton+et_neutrino),2) - pow((lepton[0]+neutrino[0]),2) - pow((lepton[1]+neutrino[1]),2)));
                                InvMassLPN->Fill(sqrt(pow((lepton[4]+photon[4]+neutrino[4]),2) - pow((lepton[0]+photon[0]+neutrino[0]),2) - pow((lepton[1]+photon[1]+neutrino[1]),2) - pow((lepton[2]+photon[2]+neutrino[2]),2)));
                                InvMassLN->Fill(sqrt(pow((lepton[4]+neutrino[4]),2) - pow((lepton[0]+neutrino[0]),2) - pow((lepton[1]+neutrino[1]),2) - pow((lepton[2]+neutrino[2]),2)));
                                MissingE->Fill(sqrt( pow((neutrino1[4]+neutrino2[4]+neutrino3[4]),2) - pow((neutrino1[0]+neutrino2[0]+neutrino3[0]),2) - pow((neutrino1[1]+neutrino2[1]+neutrino3[1]),2)));
                                MEwLepton->Fill(neutrino[4]);
                                double DeltaPhiLP_val = fabs(lepton[6]-photon[6]);
                                if(DeltaPhiLP_val > pi){
                                    DeltaPhiLP_val = 2*pi - DeltaPhiLP_val;
                                }
                                DeltaPhiLP->Fill(DeltaPhiLP_val);
                                double DeltaPhiLN_val = fabs(lepton[6]-neutrino[6]);
                                if(DeltaPhiLN_val > pi){
                                    DeltaPhiLN_val = 2*pi - DeltaPhiLN_val;
                                }
                                DeltaPhiLN->Fill(DeltaPhiLN_val);
                                DeltaEtaLP->Fill(fabs(lepton[7]-photon[7]));
                                DeltaEtaLN->Fill(fabs(lepton[7]-neutrino[7]));
                                DeltaRLP->Fill(sqrt(pow(DeltaPhiLP_val,2) + pow((lepton[7]-photon[7]),2)));
                                DeltaRLN->Fill(sqrt(pow(DeltaPhiLN_val,2) + pow((lepton[7]-neutrino[7]),2)));
                                
                                LeptonPhi->Fill(lepton[6]);
                                LeptonEta->Fill(fabs(lepton[7]));
                                PhotonPhi->Fill(photon[6]);
                                PhotonEta->Fill(fabs(photon[7]));
                                
                                //                These are the histograms for the Singly Charged and the Uncharged Higgs
                                if(photon3[3]>0){
                                    if(photon_lowest_deltar==1){
                                        double DeltaPhi2Photons_val = fabs(photon2[6]-photon3[6]);
                                        if(DeltaPhi2Photons_val > pi){
                                            DeltaPhi2Photons_val = 2*pi - DeltaPhi2Photons_val;
                                        }
                                        DeltaPhi2Photons->Fill(DeltaPhi2Photons_val);
                                        DeltaEta2Photons->Fill(fabs(photon2[7]-photon3[7]));
                                    
                                        DeltaR2Photons->Fill(sqrt(pow(DeltaPhi2Photons_val,2)+pow((photon2[7]-photon3[7]),2)));
                                        InvMass2Photons->Fill(sqrt(pow((photon2[4]+photon3[4]),2) - pow((photon2[0]+photon3[0]),2) - pow((photon2[1]+photon3[1]),2) - pow((photon2[2]+photon3[2]),2)));
                                        Mt2Photons->Fill(sqrt(pow((photon2[3]+photon3[3]),2) - pow((photon2[0]+photon3[0]),2) - pow((photon2[1]+photon3[1]),2)));
                                        Photon1E->Fill(photon1[3]);
                                        Photon2E->Fill(photon2[3]);
                                        Photon2E->Fill(photon3[3]);
                                    }else if(photon_lowest_deltar==2){
                                        double DeltaPhi2Photons_val = fabs(photon1[6]-photon3[6]);
                                        if(DeltaPhi2Photons_val > pi){
                                            DeltaPhi2Photons_val = 2*pi - DeltaPhi2Photons_val;
                                        }
                                        DeltaPhi2Photons->Fill(DeltaPhi2Photons_val);
                                        DeltaEta2Photons->Fill(fabs(photon1[7]-photon3[7]));
                                        
                                        DeltaR2Photons->Fill(sqrt(pow(DeltaPhi2Photons_val,2)+pow((photon1[7]-photon3[7]),2)));
                                        InvMass2Photons->Fill(sqrt(pow((photon1[4]+photon3[4]),2) - pow((photon1[0]+photon3[0]),2) - pow((photon1[1]+photon3[1]),2) - pow((photon1[2]+photon3[2]),2)));
                                        Mt2Photons->Fill(sqrt(pow((photon1[3]+photon3[3]),2) - pow((photon1[0]+photon3[0]),2) - pow((photon1[1]+photon3[1]),2)));
                                        Photon1E->Fill(photon2[3]);
                                        Photon2E->Fill(photon1[3]);
                                        Photon2E->Fill(photon3[3]);
                                    }else{
                                        double DeltaPhi2Photons_val = fabs(photon2[6]-photon1[6]);
                                        if(DeltaPhi2Photons_val > pi){
                                            DeltaPhi2Photons_val = 2*pi - DeltaPhi2Photons_val;
                                        }
                                        DeltaPhi2Photons->Fill(DeltaPhi2Photons_val);
                                        DeltaEta2Photons->Fill(fabs(photon2[7]-photon1[7]));
                                        
                                        DeltaR2Photons->Fill(sqrt(pow(DeltaPhi2Photons_val,2)+pow((photon2[7]-photon1[7]),2)));
                                        InvMass2Photons->Fill(sqrt(pow((photon2[4]+photon1[4]),2) - pow((photon2[0]+photon1[0]),2) - pow((photon2[1]+photon1[1]),2) - pow((photon2[2]+photon1[2]),2)));
                                        Mt2Photons->Fill(sqrt(pow((photon2[3]+photon1[3]),2) - pow((photon2[0]+photon1[0]),2) - pow((photon2[1]+photon1[1]),2)));
                                        Photon1E->Fill(photon3[3]);
                                        Photon2E->Fill(photon2[3]);
                                        Photon2E->Fill(photon1[3]);
                                    }
                                }
                                
                                //                These are the histograms for the Singly Charged and the Doubly Charged Higgs
                                if(quark1[3]>0){
                                    double DeltaPhiJets1_val = fabs(quark1[6]-quark2[6]);
                                    if(DeltaPhiJets1_val > pi){
                                        DeltaPhiJets1_val = 2*pi - DeltaPhiJets1_val;
                                    }
                                    DeltaPhiJets1->Fill(DeltaPhiJets1_val);
                                    DeltaEtaJets1->Fill(fabs(quark1[7]-quark2[7]));
                                    DeltaRJets1->Fill(sqrt(pow(DeltaPhiJets1_val,2) + pow((quark1[7]-quark2[7]),2)));
                                    InvMassJets1->Fill(sqrt(pow((quark1[4]+quark2[4]),2) - pow((quark1[0]+quark2[0]),2) - pow((quark1[1]+quark2[1]),2) - pow((quark1[2]+quark2[2]),2)));
                                    double et_quark1 = sqrt(pow(quark1[5],2)+pow(quark1[0],2)+pow(quark1[1],2));
                                    double et_quark2 = sqrt(pow(quark2[5],2)+pow(quark2[0],2)+pow(quark2[1],2));
                                    MtJets1->Fill(sqrt(pow((et_quark1+et_quark2),2) - pow((quark1[0]+quark2[0]),2) - pow((quark1[1]+quark2[1]),2)));
                                }
                                
                                if(quark3[3]>0){
                                    double DeltaPhiJets2_val = fabs(quark3[6] - quark4[6]);
                                    if(DeltaPhiJets2_val > pi){
                                        DeltaPhiJets2_val = 2*pi - DeltaPhiJets2_val;
                                    }
                                    DeltaPhiJets2->Fill(DeltaPhiJets2_val);
                                    DeltaEtaJets2->Fill(fabs(quark3[7]-quark4[7]));
                                    DeltaRJets2->Fill(sqrt(pow(DeltaPhiJets2_val,2) + pow((quark3[7]-quark4[7]),2)));
                                    InvMassJets2->Fill(sqrt(pow((quark3[4]+quark4[4]),2) - pow((quark3[0]+quark4[0]),2) - pow((quark3[1]+quark4[1]),2) - pow((quark3[2]+quark4[2]),2)));
                                    double et_quark3 = sqrt(pow(quark3[5],2)+pow(quark3[0],2)+pow(quark3[1],2));
                                    double et_quark4 = sqrt(pow(quark4[5],2)+pow(quark4[0],2)+pow(quark4[1],2));
                                    MtJets2->Fill(sqrt(pow((et_quark3+et_quark4),2) - pow((quark3[0]+quark4[0]),2) - pow((quark3[1]+quark4[1]),2)));
                                }
                                if(lepton3[3]>0){
                                    Lepton1Pt->Fill(lepton1[3]);
                                    Lepton2Pt->Fill(lepton2[3]);
                                    Lepton2Pt->Fill(lepton3[3]);
                                }
                                if(photon2[3]>0 and photon3[3]==0.0){
                                    double DeltaPhiSepPhotons_val = fabs(photon1[6] - photon2[6]);
                                    if(DeltaPhiSepPhotons_val > pi){
                                        DeltaPhiSepPhotons_val = 2*pi - DeltaPhiSepPhotons_val;
                                    }
                                    DeltaPhiSepPhotons->Fill(DeltaPhiSepPhotons_val);
                                    DeltaEtaSepPhotons->Fill(fabs(photon1[7]-photon2[7]));
                                    DeltaRSepPhotons->Fill(sqrt(pow(DeltaPhiSepPhotons_val,2) + pow((photon1[7]-photon2[7]),2)));
                                }
                                
                                if(photon_lowest_deltar==1 and photon_passed[0]){
                                    ClosestPhoton->Fill(0);
                                }else if(photon_lowest_deltar==2 and photon_passed[1]){
                                    ClosestPhoton->Fill(1);
                                }else{
                                    ClosestPhoton->Fill(2);
                                }
                                if(lepton_highest_pt==1){
                                    HighestPtLepton->Fill(0);
                                }else if(lepton_highest_pt==2){
                                    HighestPtLepton->Fill(1);
                                }else{
                                    HighestPtLepton->Fill(2);
                                }
                                if(neutrino_lowest_deltar==1){
                                    ClosestNeutrino->Fill(0);
                                }else if(neutrino_lowest_deltar==2){
                                    ClosestNeutrino->Fill(1);
                                }else{
                                    ClosestNeutrino->Fill(2);
                                }
                                
                                if(lepton_highest_pt==1 and neutrino_lowest_deltar==1){
                                    SameNeutrino->Fill(1);
                                }else if(lepton_highest_pt==2 and neutrino_lowest_deltar==2){
                                    SameNeutrino->Fill(1);
                                }else if(lepton_highest_pt==3 and neutrino_lowest_deltar==3){
                                    SameNeutrino->Fill(1);
                                }else{
                                    SameNeutrino->Fill(0);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        float efficiency = num_all_cuts*100.0/num_cuts[0];
        float fid_efficiency = num_fid_cuts[4]*100.0/num_fid_cuts[0];
        
        printf("This is the number of events that make each cut: Total number: %d\nFirst cut: %d\nSecond cut: %d\nThird cut: %d\nFourth cut: %d\nAll cuts: %d\nAnd this is the efficiency: %f\n",num_cuts[0],num_cuts[1],num_cuts[2],num_cuts[3],num_cuts[4],num_all_cuts,efficiency);
        printf("These are the fiducial values: %d\nAfter pt cut: %d\nAfter photon cut: %d\nAfter jet cut: %d\nAfter b-jet cut: %d\nand this is the efficiency: %f\n",num_fid_cuts[0],num_fid_cuts[1],num_fid_cuts[2],num_fid_cuts[3],num_fid_cuts[4],fid_efficiency);
        
        std::string output_file = line.substr(0,line.size()-5);
        output_file += "_output.root";
        TFile *outfile = new TFile(output_file.c_str(), "RECREATE");

        PtLGammaNu->Write();
        MET->Write();
        HT->Write();
        Plq->Write();
        
        //        These are general histograms that fit every signal
        ParticleCounter->Write();
        LeptonPt->Write();
        LeptonP->Write();
        PhotonEt->Write();
        PhotonE->Write();
        MtLPN->Write();
        MtLN->Write();
        InvMassLPN->Write();
        InvMassLN->Write();
        
        MissingE->Write();
        MEwLepton->Write();
        
        DeltaPhiLP->Write();
        DeltaPhiLN->Write();
        DeltaEtaLP->Write();
        DeltaEtaLN->Write();
        DeltaRLP->Write();
        DeltaRLN->Write();
        
        LeptonPhi->Write();
        LeptonEta->Write();
        PhotonPhi->Write();
        PhotonEta->Write();
        
        //        These are histograms that only apply to the first signal (charged and uncharged higgs)
        DeltaPhi2Photons->Write();
        DeltaEta2Photons->Write();
        DeltaR2Photons->Write();
        ClosestPhoton->Write();
        InvMass2Photons->Write();
        Mt2Photons->Write();
        
        //        These are histograms that only apply to the second signal (charged and doubly charged higgs)
        DeltaPhiJets1->Write();
        DeltaPhiJets2->Write();
        DeltaEtaJets1->Write();
        DeltaEtaJets2->Write();
        DeltaRJets1->Write();
        DeltaRJets2->Write();
        InvMassJets1->Write();
        InvMassJets2->Write();
        MtJets1->Write();
        MtJets2->Write();
        
        //        These are histograms that only apply to the third signal (two charged higgs)
        DeltaPhiSepPhotons->Write();
        DeltaEtaSepPhotons->Write();
        DeltaRSepPhotons->Write();
        
        Lepton1Pt->Write();
        Lepton2Pt->Write();
        
        Photon1E->Write();
        Photon2E->Write();
        
        HighestPtLepton->Write();
        ClosestNeutrino->Write();
        
        SameNeutrino->Write();
        
        outfile->Close();
        myROOTFile->Close();
        printf("This is the time it takes to analyse the particular file: %f\n",double((clock()-start_time)/1e6));
    }
//    HistOverlay();
    printf("This is the time it takes to analyse the particular file: %f\n",double((clock()-code_start_time)/1e6));
}


int maximum(float a, float b, float c){
    if(a>b){
        if(a>c){
            return 1;
        }else{
            return 3;
        }
    }else{
        if(b>c){
            return 2;
        }else{
            return 3;
        }
    }
}

int minimum(float a, float b, float c){
    if(a<b){
        if(a<c){
            return 1;
        }else{
            return 3;
        }
    }else{
        if(b<c){
            return 2;
        }else{
            return 3;
        }
    }
}

int jet_cuts(float a, float b, float c, float d, float e, float f){
    int num_jets = 0;
    if(a>20.0){
        num_jets += 1;
    }
    if(b>20.0){
        num_jets += 1;
    }
    if(c>20.0){
        num_jets += 1;
    }
    if(d>20.0){
        num_jets += 1;
    }
    if(e>20.0){
        num_jets += 1;
    }
    if(f>20.0){
        num_jets += 1;
    }
    return num_jets;
}
