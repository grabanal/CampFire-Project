//
//  HistOverlay.cpp
//
//
//  Created by Utsav Patel on 5/6/19.
//

#include <stdio.h>
#include "TRandom3.h"
#include "TChain.h"
#include "TProfile.h"
#include "TCut.h"
#include "TLegend.h"
#include <cmath>
#include <cstdlib>
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TMemberInspector.h"
#include "TObjArray.h"
#include "TSeqCollection.h"
#include "TNamed.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom.h>
#include <vector>
#include <exception>
#include <string.h>
#include "TLatex.h"
#include "TDirectory.h"
#include <algorithm>
#include "TPad.h"
#include "TLegend.h"

void HistogramOverlay(){
    
    TFile infile1("BG_output.root");
    TH1F *h1_MET = (TH1F*)infile1.Get("MET");
    TH1F *h1_HT = (TH1F*)infile1.Get("HT");
    TH1F *h1_PtLGammaNu = (TH1F*)infile1.Get("PtLGammaNu");
    TH1F *h1_Plq = (TH1F*)infile1.Get("Plq");
    
    TFile infile2("Signal_output.root");
    TH1F *h2_MET = (TH1F*)infile2.Get("MET");
    TH1F *h2_HT = (TH1F*)infile2.Get("HT");
    TH1F *h2_PtLGammaNu = (TH1F*)infile2.Get("PtLGammaNu");
    TH1F *h2_Plq = (TH1F*)infile2.Get("Plq");
    
    double h1_MET_Scale = 1./(h1_MET->Integral());
    double h2_MET_Scale = 1./(h2_MET->Integral());
    
    double h1_HT_Scale = 1./(h1_HT->Integral());
    double h2_HT_Scale = 1./(h2_HT->Integral());
    
    double h1_PtLGammaNu_Scale = 1./(h1_PtLGammaNu->Integral());
    double h2_PtLGammaNu_Scale = 1./(h2_PtLGammaNu->Integral());
    
    double h1_Plq_Scale = 1./(h1_Plq->Integral());
    double h2_Plq_Scale = 1./(h2_Plq->Integral());
    
    h1_MET->Scale(h1_MET_Scale);
    h2_MET->Scale(h2_MET_Scale);
    
    h1_HT->Scale(h1_HT_Scale);
    h2_HT->Scale(h2_HT_Scale);
    
    h1_PtLGammaNu->Scale(h1_PtLGammaNu_Scale);
    h2_PtLGammaNu->Scale(h2_PtLGammaNu_Scale);
    
    h1_Plq->Scale(h1_Plq_Scale);
    h2_Plq->Scale(h2_Plq_Scale);
    
    TCanvas *c1 = new TCanvas("c1", "Missing E_{T}", 600, 400);
    TCanvas *c2 = new TCanvas("c2", "H_{T}", 600, 400);
    TCanvas *c3 = new TCanvas("c3", "P_{T}^{l+#gamma+MET}", 600, 400);
    TCanvas *c4 = new TCanvas("c4", "P_{l} #dot q", 600, 400);
    
    c1->cd(0);
    h1_MET->Draw("HIST same");
    h1_MET->SetLineColor(kRed);
    h2_MET->Draw("HIST same");
    h2_MET->SetLineColor(kBlue);
    auto* legend_1 = new TLegend(0.76,0.70,0.98,0.97);
    TLegendEntry *h1_MET_entry = legend_1->AddEntry("h1_MET","W#gamma","l");
    h1_MET_entry->SetLineColor(kRed);
    TLegendEntry *h2_MET_entry = legend_1->AddEntry("h2_MET","H^{#pm}H^{#mp#mp}","l");
    h2_MET_entry->SetLineColor(kBlue);
    legend_1->Draw("same");
    gPad->Modified();
    gPad->Update();
    c1->Update();
    c1->SaveAs("MET.pdf");
    
    c2->cd(0);
    h1_HT->Draw("HIST same");
    h1_HT->SetLineColor(kRed);
    h2_HT->Draw("HIST same");
    h2_HT->SetLineColor(kBlue);
    auto* legend_2 = new TLegend(0.76,0.70,0.98,0.97);
    TLegendEntry *h1_HT_entry = legend_2->AddEntry("h1_HT","W#gamma","l");
    h1_HT_entry->SetLineColor(kRed);
    TLegendEntry *h2_HT_entry = legend_2->AddEntry("h2_HT","H^{#pm}H^{#mp#mp}","l");
    h2_HT_entry->SetLineColor(kBlue);
    legend_2->Draw("same");
    gPad->Modified();
    gPad->Update();
    c2->Update();
    c2->SaveAs("HT.pdf");
    
    c3->cd(0);
    h1_PtLGammaNu->Draw("HIST same");
    h1_PtLGammaNu->SetLineColor(kRed);
    h2_PtLGammaNu->Draw("HIST same");
    h2_PtLGammaNu->SetLineColor(kBlue);
    auto* legend_3 = new TLegend(0.76,0.70,0.98,0.97);
    TLegendEntry *h1_PtLGammaNu_entry = legend_3->AddEntry("h1_PtLGammaNu","W#gamma","l");
    h1_PtLGammaNu_entry->SetLineColor(kRed);
    TLegendEntry *h2_PtLGammaNu_entry = legend_3->AddEntry("h2_PtLGammaNu","H^{#pm}H^{#mp#mp}","l");
    h2_PtLGammaNu_entry->SetLineColor(kBlue);
    legend_3->Draw("same");
    gPad->Modified();
    gPad->Update();
    c3->Update();
    c3->SaveAs("PTLGammaNu.pdf");
    
    c4->cd(0);
    h1_Plq->Draw("HIST same");
    h1_Plq->SetLineColor(kRed);
    h2_Plq->Draw("HIST same");
    h2_Plq->SetLineColor(kBlue);
    auto* legend_4 = new TLegend(0.76,0.70,0.98,0.97);
    TLegendEntry *h1_Plq_entry = legend_4->AddEntry("h1_Plq","W#gamma","l");
    h1_Plq_entry->SetLineColor(kRed);
    TLegendEntry *h2_Plq_entry = legend_4->AddEntry("h2_Plq","H^{#pm}H^{#mp#mp}","l");
    h2_Plq_entry->SetLineColor(kBlue);
    legend_4->Draw("same");
    gPad->Modified();
    gPad->Update();
    c4->Update();
    c4->SaveAs("Plq.pdf");
}

