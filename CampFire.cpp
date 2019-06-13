//
//  CampFire.cpp
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
int jet_cuts(float,float,float,float);
//int * photon_filter(int,int,int);

void CampFire(){
    
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
        TRandom *rand_num = new TRandom();
        int num_events = particle_Tree->GetEntries();
        int num_cuts[4] = {0};
        
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
            
            float lepton[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //This is the lepton with the highest pt
            float photon[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //This is the closest photon to the lepton with the highest pt
            
            bool lepton_present[2] = {false,false};
            bool neutrino_present[3] = {false,false,false};
            bool quark_present[3] = {false,false,false};
            bool photon_passed = false;
            
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
                    photon[0] = px->GetValue(j);
                    photon[1] = py->GetValue(j);
                    photon[2] = pz->GetValue(j);
                    photon[3] = pt->GetValue(j);
                    photon[4] = energy->GetValue(j);
                    photon[5] = mass->GetValue(j);
                    photon[6] = phi->GetValue(j);
                    photon[7] = eta->GetValue(j);
                    if(photon[3]>25.0 and fabs(photon[7])<2.5){
                        photon_passed = true;
                    }
                }else if((particleID==1 || particleID==2 || particleID==3 || particleID==4 || particleID==-1 || particleID==-2 || particleID==-3 || particleID==-4) and j>1){
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
            
            num_cuts[0] += 1;
            
            float missing_et = sqrt(pow((neutrino1[0]+neutrino2[0]+neutrino3[0]),2) + pow((neutrino1[1]+neutrino2[1]+neutrino3[1]),2)); //Assume neutrino mass is approximately zero
            float h_t = lepton1[3]+lepton2[3]+lepton3[3]+quark1[3]+quark2[3]+quark3[3]+quark4[3]+photon[3];
            float ptlgnu = sqrt(pow((lepton[0]+photon[0]),2) + pow((lepton[1]+photon[1]),2) + pow((neutrino1[0]+neutrino2[0]+neutrino3[0]),2) + pow((neutrino1[1]+neutrino2[1]+neutrino3[1]),2));
            float pl_q = lepton[4]*photon[4] - (lepton[0]*photon[0]+lepton[1]*photon[1]+lepton[2]*photon[2]);
            
            //            First cut
            if(lepton[3]>25.0 and fabs(lepton[7])<2.5){
                num_cuts[1] += 1;
                made_cuts[0] = true;
            }
            //            Second cut
            if(photon_passed){
                num_cuts[2] += 1;
                made_cuts[1] = true;
            }
            //            Third cut
            int num_qualifying_jets = jet_cuts(quark1[3],quark2[3],quark3[3],quark4[3]);
            if(num_qualifying_jets<=2){
                num_cuts[3] += 1;
                made_cuts[2] = true;
            }
            
            if((made_cuts[0] && made_cuts[1]) && (made_cuts[2])){
                num_all_cuts += 1;
                //                These are where I add every other plot. These can be moved to the fiducial cuts block if need be
                PtLGammaNu->Fill(ptlgnu);
                MET->Fill(missing_et);
                HT->Fill(h_t);
                Plq->Fill(pl_q);
            }
        }
        
        float efficiency = num_all_cuts*100.0/num_cuts[0];
        
        printf("This is the number of events that make each cut: Total number: %d\nFirst cut: %d\nSecond cut: %d\nThird cut: %d\nAll cuts: %d\nAnd this is the efficiency: %f\n",num_cuts[0],num_cuts[1],num_cuts[2],num_cuts[3],num_all_cuts,efficiency);
        
        std::string output_file = line.substr(0,line.size()-5);
        output_file += "_output.root";
        TFile *outfile = new TFile(output_file.c_str(), "RECREATE");
        
        PtLGammaNu->Write();
        MET->Write();
        HT->Write();
        Plq->Write();
        
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

int jet_cuts(float a, float b, float c, float d){
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
    return num_jets;
}
