
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"

void eNuE(){
    TFile f("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_intrinsic_nue_overlay.root");
    
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
    t-> SetBranchAddress("trueSimPartProcess", trueSimPartProcess);

    for(Int_t i = 0; i < t->GetEntries(); ++i){

        t -> GetEntry(i);
        
        if(trueNuPDG == 12 || trueNuPDG == -12){
            for(Int_t j = 0; j < nTruePrimParts; ++j){
                if(truePrimPartPDG[j] == 211 || truePrimPartPDG[j] == -211){
                    h1.Fill(truePrimPartE[j], xsecWeight);
                    
                    for(int l = 0; l < nTrueSimParts; ++l){
                        if(trueSimPartProcess[l] == 0){
                            if(trueSimPartPDG[l] == 211 || trueSimPartPDG[l] == -211){
                                for(int k = 0; k < nTracks; ++k){
                                    if(trackTrueTID[k] == trueSimPartTID[l]){
                                        h2.Fill(truePrimPartE[j], xsecWeight);
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
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
