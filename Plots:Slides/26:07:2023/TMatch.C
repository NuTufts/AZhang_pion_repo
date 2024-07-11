
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

//fill and plot a histogram of reconstructed track energies

void TMatch(){

//  TH1F h("h","TrackTrueComp for Confirmed Pions", 40,0,1); //create a histogram for reconstructed track energies

  TFile f("prepare_selection_test_reco_v2me05_gen2val_v23_nu_file.root");
  TTree* t = (TTree*)f.Get("EventTree");

  Int_t nTracks; //the nTracks values will get assigned to this variable when looping through the tree. it is an int
  t -> SetBranchAddress("nTracks", &nTracks); // first arg: variable name, second arg: memory address for variable
//  Float_t trackTrueComp[100]; //the track reconstructed energy values will get assigned to this variable when looping through the tree. it is an array of floating point variables, with one entry for each track in the event. Choose a size for the array here (100) that will be larger than the actual array size for any particular event to make sure it is large enough to always accomodate all track values
//  t -> SetBranchAddress("trackTrueComp", trackTrueComp); //since this variable is for an array, you don't need the "&" in the second argument, array variables just point to the memory address in C++
//
  Int_t trackPID[100];
  Int_t trackTruePID[100];
  t -> SetBranchAddress("trackPID", trackPID);
  t -> SetBranchAddress("trackTruePID", trackTruePID);
    
    int pi_match, electron_match, muon_match, proton_match, csm_ray, t_total = 0;

  for(int i = 0; i < t->GetEntries(); ++i){

    t -> GetEntry(i); //this will assign the number of tracks for entry i to the nTracks variable we defined above and similarly for the trackRecoE array

      for(int j = 0; j < nTracks; ++j){
          if(trackTruePID[j] == 211){
              t_total += 1;
          }
          if(trackPID[j] == 211 && trackTruePID[j] == 211){
              pi_match += 1;
          }
          if(trackPID[j] == 11 && trackTruePID[j] == 211){
              electron_match += 1;
          }
          if(trackPID[j] == 13 && trackTruePID[j] == 211){
              muon_match += 1;
          }
          if(trackPID[j] == 2212 && trackTruePID[j] == 211){
              proton_match += 1;
              //          }else
              //          if(trackPID[j] == 22 && trackTruePID[j] == 211){
              //              gamma += 1;
          }
          if(trackPID[j] == 0 && trackTruePID[j] == 211){
              csm_ray += 1;
          }
      }
    
  }
    cout << "pion matches: " << pi_match << endl
         << "electron matches: " << electron_match << endl
         << "muon matches: " << muon_match << endl
         << "proton matches: " << proton_match << endl
         << "cosmic ray matches: " << csm_ray << endl
         << "total tracks: " << t_total << endl;
//  
//  h.Draw("EHIST"); //this should make a plot of the histogram pop up on your screen, but it doesn't seem to be working for me, so write it to a new output file as well:
//  TFile f_out("TrueComp.root","RECREATE");
//  h.Write();
//  //if you open this output_plots.root file in a TBrowser after running the script, you can double click on the histogram to plot it

}
