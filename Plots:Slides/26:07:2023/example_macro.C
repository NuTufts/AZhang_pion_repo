
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

//fill and plot a histogram of reconstructed track energies

void example_macro(){

  TH1F h("h","reconstructed track energy", 100,0,1500); //create a histogram for reconstructed track energies

  TFile f("prepare_selection_test_reco_v2me05_gen2val_v23_nu_file.root");
  TTree* t = (TTree*)f.Get("EventTree");

  Int_t nTracks; //the nTracks values will get assigned to this variable when looping through the tree. it is an int
  t -> SetBranchAddress("nTracks", &nTracks); // first arg: variable name, second arg: memory address for variable
  Float_t trackRecoE[100]; //the track reconstructed energy values will get assigned to this variable when looping through the tree. it is an array of floating point variables, with one entry for each track in the event. Choose a size for the array here (100) that will be larger than the actual array size for any particular event to make sure it is large enough to always accomodate all track values
  t -> SetBranchAddress("trackRecoE", trackRecoE); //since this variable is for an array, you don't need the "&" in the second argument, array variables just point to the memory address in C++

  for(Int_t i = 0; i < t->GetEntries(); ++i){

    t -> GetEntry(i); //this will assign the number of tracks for entry i to the nTracks variable we defined above and similarly for the trackRecoE array

    for(Int_t j = 0; j < nTracks; ++j){
      h.Fill(trackRecoE[j]); //add the reconstructed energy of the jth track in this event to the histogram
    }

  }
  
  h.Draw("EHIST"); //this should make a plot of the histogram pop up on your screen, but it doesn't seem to be working for me, so write it to a new output file as well:
  TFile f_out("output_plots.root","RECREATE");
  h.Write();
  //if you open this output_plots.root file in a TBrowser after running the script, you can double click on the histogram to plot it

}
