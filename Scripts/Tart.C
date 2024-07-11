
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"

void eNuE(){
    TFile f("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_nu_overlay.root");
    
    int t_pot = 0;
    Float_t totGoodPOT;
    TTree* p = (TTree*)f.Get("potTree");
    p -> SetBranchAddress("totGoodPOT", &totGoodPOT);
    for(Int_t i = 0; i < p->GetEntries(); ++i){
        p->GetEntry(i);
        t_pot += totGoodPOT;
    }
    
    TH1F h1("h1","Simulated Pion Energy from Electron Neutrino Interactions;TruePrimPartE (GeV);# of Pions", 40,0,4);

    TH1F h2("h2","Truth-Matched Pion Energy from Electron Neutrino Interactions;TruePrimPartE (GeV);# of Pions", 40,0,4);

    TTree* t = (TTree*)f.Get("EventTree");
    Float_t xsecWeight;
    Int_t trueNuPDG;
    Int_t nTruePrimParts;
    Int_t nTracks;
    Int_t trackTrueTID[200];
    Int_t nTrueSimParts;
    Int_t trueSimPartTID[200];
    Float_t truePrimPartE[200];
    Int_t truePrimPartPDG[200];
    Int_t trueSimPartPDG[200];
    Int_t trueSimPartProcess[200];
    Float_t trueNuE;
    
    t-> SetBranchAddress("trueSimPartPDG", trueSimPartPDG);
    t-> SetBranchAddress("trueNuPDG", &trueNuPDG);
    t-> SetBranchAddress("nTruePrimParts", &nTruePrimParts);
    t-> SetBranchAddress("nTracks", &nTracks);
    t-> SetBranchAddress("trackTrueTID", trackTrueTID);
    t-> SetBranchAddress("nTrueSimParts", &nTrueSimParts);
    t-> SetBranchAddress("trueSimPartTID", trueSimPartTID);
    t-> SetBranchAddress("truePrimPartE", truePrimPartE);
    t-> SetBranchAddress("truePrimPartPDG", truePrimPartPDG);
    t-> SetBranchAddress("xsecWeight", &xsecWeight);
    t-> SetBranchAddress("trueNuE", &trueNuE);
    t-> SetBranchAddress("trueSimPartProcess", trueSimPartProcess);
    
    for(Int_t i = 0; i < t->GetEntries(); ++i){
        
        t -> GetEntry(i);
        int pi_num = 0, mu_num = 0;
        
        
        for(Int_t j = 0; j < nTrueSimParts; ++j){
            if(trueSimPartPDG[j] == 211 || trueSimPartPDG[j] == -211){
                pi_num += 1;
            }
            if(trueSimPartPDG[j] == 13 || trueSimPartPDG[j] == -13){
                mu_num += 1;
            }
        }
        if(pi_num == 1 && mu_num == 1){
            h1.Fill(trueNuE);
        }
    }
    TEfficiency h3(h2, h1);
    h3.SetUseWeightedEvents(true);
    h3.SetTitle("Reconstruction Efficiency of Pions in Electron Neutrino Events by Energy");
    
    h1.Scale(6.67e20/t_pot);
    h2.Scale(6.67e20/t_pot);
    h1.Draw("EHIST");
    h2.Draw("EHIST");
    h3.Draw("EHIST");
    
    TFile f_out("eNuE.root","RECREATE");
    h1.Write();
    h2.Write();
    h3.Write();
}
