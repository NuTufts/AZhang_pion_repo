
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"

//fill and plot a histogram of reconstructed track energies

void pimu(){
    TFile f("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_nu_overlay.root");
    
    int t_pot = 0;
    Float_t totGoodPOT;
    TTree* p = (TTree*)f.Get("potTree");
    p -> SetBranchAddress("totGoodPOT", &totGoodPOT);
    for(Int_t i = 0; i < p->GetEntries(); ++i){
        p->GetEntry(i);
        t_pot += totGoodPOT;
    }
//    TH1F h1("h1","Simulated Pion Energy from Electron Interactions;TruePrimPartE;# of Pions", 40,0,4); //create a histogram for reconstructed track energies
//
//    TH1F h2("h2","Truth-Matched Pion Energy from Electron Interactions;TruePrimPartE;# of Pions", 40,0,4); //create a histogram for reconstructed track energies

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
    Float_t trackTrueE[200];
    Int_t trueSimPartProcess[200];
    Int_t trackTruePID[200];
    Int_t trackPID[200];
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
    t-> SetBranchAddress("trackTrueE", trackTrueE);
    t-> SetBranchAddress("trueSimPartProcess", trueSimPartProcess);
    t-> SetBranchAddress("trackTruePID", trackTruePID);
    t-> SetBranchAddress("trackPID", trackPID);

    int total = 0;
    int passed = 0;
    bool pi, mu;
    for(Int_t i = 0; i < t->GetEntries(); ++i){

        t -> GetEntry(i);
        
        for(Int_t j = 0; j < nTruePrimParts; ++j){
            if(truePrimPartPDG[j] == 211 || truePrimPartPDG[j] == -211){
                pi = true;
            }
            if(truePrimPartPDG[j] == 13 || truePrimPartPDG[j] == -13){
                mu = true;
            }
        }
        if(pi == true && mu == true){
            total += 1;
            pi = false;
            mu = false;
        }
        for(int l = 0; l < nTracks; ++l){
            if(trackTruePID[l] == trackPID[l]){
                if(trackPID[l] == 211 || trackPID[l] == -211){
                    pi = true;
                }
                if(trackPID[l] == 13 || trackPID[l] == -13){
                    mu = true;
                }
            }
        }
        if(pi == true && mu == true){
            passed += 1;
            pi = false;
            mu = false;
        }
    }
    
    cout << total << endl;
    cout << passed << endl;
    
//    TEfficiency h3(h2, h1);
////    h3.SetUseWeightedEvents(true);
//    h1.Scale(6.67e20/t_pot);
//    h2.Scale(6.67e20/t_pot);
//    h1.Draw("EHIST");
//    h2.Draw("EHIST");
////    TH1F h3 = (h2/h1);
//    TFile f_out("eNuE.root","RECREATE");
//    h1.Write();
//    h2.Write();
//    h3.Write();
    //if you open this output_plots.root file in a TBrowser after running the script, you can double click on the histogram to plot it

}
