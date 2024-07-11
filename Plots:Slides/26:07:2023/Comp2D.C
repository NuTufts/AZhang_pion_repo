
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TClass.h"
//fill and plot a histogram of reconstructed track energies

void Comp2D(){

  TH2F h("h","Comp vs. True Comp of Truth-Matched Pions;trackTrueComp (%);trackComp (%)", 40,0,1,40,0,1); //create a histogram for reconstructed track energies
    
  TFile f("prepare_selection_test_reco_v2me05_gen2val_v23_nu_file.root");
  TTree* t = (TTree*)f.Get("EventTree");

  Int_t nTracks; //the nTracks values will get assigned to this variable when looping through the tree. it is an int
  t -> SetBranchAddress("nTracks", &nTracks); // first arg: variable name, second arg: memory address for variable
  Float_t trackTrueComp[100]; //the track reconstructed energy values will get assigned to this variable when looping through the tree. it is an array of floating point variables, with one entry for each track in the event. Choose a size for the array here (100) that will be larger than the actual array size for any particular event to make sure it is large enough to always accomodate all track values
  t -> SetBranchAddress("trackTrueComp", trackTrueComp); //since this variable is for an array, you don't need the "&" in the second argument, array variables just point to the memory address in C++
  Float_t trackComp[100];
  t -> SetBranchAddress("trackComp", trackComp);


  Int_t trackPID[100];
  Int_t trackTruePID[100];
  t -> SetBranchAddress("trackPID", trackPID);
  t -> SetBranchAddress("trackTruePID", trackTruePID);

  for(Int_t i = 0; i < t->GetEntries(); ++i){

    t -> GetEntry(i); //this will assign the number of tracks for entry i to the nTracks variable we defined above and similarly for the trackRecoE array

    for(Int_t j = 0; j < nTracks; ++j){
        if(trackPID[j] == 211 || trackTruePID[j] == 211){
            h.Fill(trackTrueComp[j], trackComp[j]);
        }
    }
  }
    
//    auto c=new TCanvas("Canvas","Canvas",800,800); c->Divide(2,2);
//    c->cd(1); h.Draw("Cont1");
//    c->cd(2); h.Draw("Colz");
//    c->cd(3); h.Draw("Lego2");
//    c->cd(4); h.Draw("Surf3");
//       // Profiles and Projections
//    auto c2=new TCanvas("Canvas2","Canvas2",800,800); c2->Divide(2,2);
//    c2->cd(1); h.ProjectionX()->Draw();
//    c2->cd(2); h.ProjectionY()->Draw();
//    c2->cd(3); h.ProfileX()->Draw();
//    c2->cd(4); h.ProfileY()->Draw();
  
  h.Draw("Colz"); //this should make a plot of the histogram pop up on your screen, but it doesn't seem to be working for me, so write it to a new output file as well:
  TFile f_out("Comp2D.root","RECREATE");
  h.Write();
  //if you open this output_plots.root file in a TBrowser after running the script, you can double click on the histogram to plot it

}
