
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include <vector>

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <math.h>
// ROOT
#include "TH2F.h"
#include "TStyle.h"
// larutil
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/ClockConstants.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
// larcv
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventPGraph.h"

// larlite
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mctruth.h"



// main fuction to load in a dlana root file, make event display image, get some truth.
int main(int nargs, char** argv){
    struct pion_container{
        int run, subrun, event = -1;
        int pi_len = 0;
        bool reco_found = false;
        float top_recoComp = -1, true_pionE = -1;
        
        void toggle_reco_found(){
            if(reco_found == true){
                reco_found = false;
            }else reco_found = true;
        }
        void set_top_recoComp(int recoComp){
            top_recoComp = recoComp;
        }
        void set_true_pionE(int pionE){
            true_pionE = pionE;
        }
    };
//    std::cout << "Hello world " << "\n";
//    std::cout << "Args Required, Order Matters:\n";
//    std::cout << "Inputfile with Truth tree\n";

    int start_entry = 0;
    
    if (nargs < 3){
        std::cout << "Not Enough Args\n";
        return 1;
    }
    std::string input_file = argv[1];
    TFile f(argv[2]);
    t = f.Get("KPSRecoManagerTree");
    
    for (int entry = start_entry; entry < f.get_entries(); entry++){
        if 
    }
    

    //createlarlite storage Manager
    larlite::storage_manager* io_larlite  = new larlite::storage_manager(larlite::storage_manager::kREAD);
    // load in and intialize larlite products
    io_larlite->add_in_filename(input_file);
    io_larlite->open();

    int nentries_mc_cv = io_larlite->get_entries();

//    std::cout << "Entries in File:  " << nentries_mc_cv << "\n";

    // loop through all the entries in the file
      for (int entry=start_entry; entry < nentries_mc_cv; entry++){

          io_larlite->go_to(entry);

      //print out r,s,e
      int _run    = io_larlite->run_id();
      int _subrun = io_larlite->subrun_id();
      int _event  = io_larlite->event_id();
//      std::cout<<"(r,s,e)="<<_run<<","<<_subrun<<","<<_event<<std::endl;

      //load in input tree
      larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)io_larlite->get_data(larlite::data::kMCTruth,  "generator" );

      // get event info

     // cc or nc event?
     bool ccevent = false;
     int ccnc = ev_mctruth->at(0).GetNeutrino().CCNC();
     if (ccnc ==0 ) ccevent=true;
     if (ccevent) std::cout<<"Is a CC Event"<<std::endl;
     else std::cout<<"Is a NC Event"<<std::endl;
          
  // type of Neutrino
  int nu_pdg = ev_mctruth->at(0).GetNeutrino().Nu().PdgCode();
  if (nu_pdg== 12) std::cout<<"Muon Neutrino event "<<std::endl;
  else if (nu_pdg== -12) std::cout<<"Muon Anti Neutrino event "<<std::endl;
  else if (nu_pdg== 14) std::cout<<"Electron Neutrino event "<<std::endl;
  else if (nu_pdg== -14) std::cout<<"Electon Anti Neutrino event "<<std::endl;

  int num_protons = 0;
  int num_neutrons = 0;
  int num_pion_charged = 0;
  int num_pion_neutral = 0;
  for(int part =0;part<(int)ev_mctruth->at(0).GetParticles().size();part++){
    // pick only final state particles
    if (ev_mctruth->at(0).GetParticles().at(part).StatusCode() == 1){
      int pdg = ev_mctruth->at(0).GetParticles().at(part).PdgCode();
      if (pdg == 2212) num_protons++;
      else if (pdg == 2112) num_neutrons++;
      else if (pdg == 111 ) num_pion_neutral++;
      else if (pdg == 211 || pdg == -211) num_pion_charged++;


