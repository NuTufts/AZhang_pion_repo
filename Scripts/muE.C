
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"

//fill and plot a histogram of reconstructed track energies

void muE(){
    TFile f("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_intrinsic_nue_overlay.root");
    
    int t_pot = 0;
    Float_t totGoodPOT;
    TTree* p = (TTree*)f.Get("potTree");
    p -> SetBranchAddress("totGoodPOT", &totGoodPOT);
    for(Int_t i = 0; i < p->GetEntries(); ++i){
        p->GetEntry(i);
        t_pot += totGoodPOT;
    }
    TH1F h1("h1","Simulated Pion Energy from Muon Interactions;TruePrimPartE;# of Muons", 40,0,4); //create a histogram for reconstructed track energies

    TH1F h2("h2","Truth-Matched Pion Energy from Muon Interactions;TruePrimPartE;# of Muons", 40,0,4); //create a histogram for reconstructed track energies

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

    for(Int_t i = 0; i < t->GetEntries(); ++i){

        t -> GetEntry(i);
        
        if(trueNuPDG == 13 || trueNuPDG == -13){
            for(Int_t j = 0; j < nTruePrimParts; ++j){
                if(truePrimPartPDG[j] == 211 || truePrimPartPDG[j] == -211){
                    h1.Fill(truePrimPartE[j]);
                    
                    for(int l = 0; l < nTrueSimParts; ++l){
                        if(trueSimPartProcess[l] == 0){
                            if(trueSimPartPDG[l] == 211 || trueSimPartPDG[l] == -211){
                                for(int k = 0; k < nTracks; ++k){
                                    if(trackTrueTID[k] == trueSimPartTID[l]){
                                        h2.Fill(truePrimPartE[j]);
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
//    h3->GetYaxis->SetTitle("Reconstruction Efficiency (%)");
    h3.SetTitle("Reconstruction Efficiency of Pions in Muon Events by Energy");
    h3.Draw("EHIST");
//    h3.SetUseWeightedEvents(true);
    h1.Scale(6.67e20/t_pot);
    h2.Scale(6.67e20/t_pot);
    h1.Draw("EHIST");
    h2.Draw("EHIST");
//    TH1F h3 = (h2/h1);
    TFile f_out("eNuE.root","RECREATE");
    h1.Write();
    h2.Write();
    h3.Write();
    //if you open this output_plots.root file in a TBrowser after running the script, you can double click on the histogram to plot it

}
