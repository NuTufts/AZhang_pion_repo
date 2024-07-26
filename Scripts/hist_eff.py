# hist_eff.py
# Author: Andy Zhang
# Purpose: Find Reconstruction Efficiencies of True Primary Pions ordered by pion history,
# reco status type, and most colinear primary
# Data Created: June 5th, 2024
# Last Modified: 3rd, 2024

import ROOT as rt
import numpy as np
from ROOT import gROOT, AddressOf
f1 = rt.TFile("pion_stats.root")
f_out = rt.TFile("hist_eff.root","RECREATE")
t1 = f1.Get("pion_tree")
p1 = f1.Get("potTree")

# Count total good POT for scaling
t_pot = 0.
for entry in p1:
    t_pot += entry.totGoodPOT

# Initialize all histograms

# Reco of pions by Collinearity (cc)
h_Reco_colinear_cc = rt.TH1F("h_Reco_colinear_cc","Number of Reconstructed Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_misID_colinear_cc = rt.TH1F("h_misID_colinear_cc","Number of mis-ID'd Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Sec_colinear_cc = rt.TH1F("h_Sec_colinear_cc","Number Pions Reco'd as Secondary by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Prim_colinear_cc = rt.TH1F("h_Prim_colinear_cc","Number Pions Reco'd as Primary by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Missing_colinear_cc = rt.TH1F("h_Missing_colinear_cc","Number of Missing Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Truecolinear_cc = rt.TH1F("h_Truecolinear_cc","Number of True Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Effcolinear_cc = rt.TH1F("h_Effcolinear_cc","Reco Efficiency of Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);Efficiency (%)", 40,-1,1)
h_Stack_colinear_cc = rt.THStack("h_Stack_colinear_cc", "Reco Status of Pions by Most Colinear Primary (cc)")

# History of pions by Reco status (cc)
h_Hist_reco_cc = rt.TH1F("h_Hist_reco_cc", "History of Reconstructed Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_misID_cc = rt.TH1F("h_Hist_misID_cc", "History of misID'd Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_sec_cc = rt.TH1F("h_Hist_sec_cc", "History of Pions ID'd as Secondary (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_prim_cc = rt.TH1F("h_Hist_prim_cc", "History of Pions ID'd as Primary (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_missing_cc = rt.TH1F("h_Hist_missing_cc", "History of Missing Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
# h_Hist_stack_cc = rt.THStack("h_Hist_stack_cc", "History of Pions by Reco Type (cc)")

# Reco of pions by History (cc)
h_Uncontained_cc = rt.TH1F("h_Uncontained_cc", "Reco Result of Uncontained Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Range_cc = rt.TH1F("h_Range_cc", "Reco Result of Ranged out Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Flight_decay_cc = rt.TH1F("h_Flight_decay_cc", "Reco Result of Decayed in Flight Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Sec_int_cc = rt.TH1F("h_Sec_int_cc", "Reco Result of Secondary Interacted Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)

# Distance to Interaction of pions by Reco Status (cc)
h_Dist_reco_cc = rt.TH1F("h_Dist_reco_cc", "Distance to Interaction of Reconstructed Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_misID_cc = rt.TH1F("h_Dist_misID_cc", "Distance to Interaction of mis-ID Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_sec_cc = rt.TH1F("h_Dist_sec_cc", "Distance to Interaction of ID'd as Sec. Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_prim_cc = rt.TH1F("h_Dist_prim_cc", "Distance to Interaction of ID'd as Prim. Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_missing_cc = rt.TH1F("h_Dist_missing_cc", "Distance to Interaction of Missing Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_stack_cc = rt.THStack("h_Dist_stack_cc", "Distance to Interaction of Pions by Reco Type (cc)")
h_Dist_eff_cc = rt.TH1F("h_Dist_eff_cc","Efficiency of Pion Reco by Distance to Interaction (cc);Distance to Interaction (cm);Efficiency(%)",100,0,50)
h_Dist_true_cc = rt.TH1F("h_Dist_true_cc", "Distance to Interaction of True Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_KE_cc = rt.TH2F("h_Dist_KE_cc", "Distance to Interaction vs True Kinetic Energy of Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)

# Distance to Interaction of Pions by History and Reco status (cc)
h_Dist_reco_hist_cc = rt.THStack("h_Dist_reco_hist_cc", "Distacce to Interaction of Reconstructed Pions by History (cc)")
h_Dist_misID_hist_cc = rt.THStack("h_Dist_misID_hist_cc", "Distacce to Interaction of mis-ID Pions by History (cc)")
h_Dist_sec_hist_cc = rt.THStack("h_Dist_sec_hist_cc", "Distacce to Interaction of ID'd as Sec. Pions by History (cc)")
h_Dist_prim_hist_cc = rt.THStack("h_Dist_prim_hist_cc", "Distance to Interaction of ID'd as Prim. Pions by History (cc)")
h_Dist_missing_hist_cc = rt.THStack("h_Dist_missing_hist_cc", "Distacce to Interaction of Missing Pions by History (cc)")
Dist_hist_reco_list_cc = [h_Dist_reco_hist_cc,h_Dist_misID_hist_cc,h_Dist_sec_hist_cc,h_Dist_prim_hist_cc,h_Dist_missing_hist_cc]

# Reco of pions by Collinearity (cc)
h_Reco_colinear_nc = rt.TH1F("h_Reco_colinear_nc","Number of Reconstructed Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_misID_colinear_nc = rt.TH1F("h_misID_colinear_nc","Number of mis-ID'd Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Sec_colinear_nc = rt.TH1F("h_Sec_colinear_nc","Number Pions Reco'd as Secondary by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Prim_colinear_nc = rt.TH1F("h_Prim_colinear_nc","Number Pions Reco'd as Primary by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Missing_colinear_nc = rt.TH1F("h_Missing_colinear_nc","Number of Missing Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Truecolinear_nc = rt.TH1F("h_Truecolinear_nc","Number of True Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Effcolinear_nc = rt.TH1F("h_Effcolinear_nc","Reco Efficiency of Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);Efficiency (%)", 40,-1,1)
h_Stack_colinear_nc = rt.THStack("h_Stack_colinear_nc", "Reco Status of Pions by Most Colinear Primary (nc)")

# History of pions by Reco status (cc)
h_Hist_reco_nc = rt.TH1F("h_Hist_reco_nc", "History of Reconstructed Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_misID_nc = rt.TH1F("h_Hist_misID_nc", "History of misID'd Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_sec_nc = rt.TH1F("h_Hist_sec_nc", "History of Pions ID'd as Secondary (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_prim_nc = rt.TH1F("h_Hist_prim_nc", "History of Pions ID'd as Primary (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_missing_nc = rt.TH1F("h_Hist_missing_nc", "History of Missing Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
# h_Hist_stack_nc = rt.THStack("h_Hist_stack_nc", "History of Pions by Reco Type (cc)")

# Reco of pions by History (cc)
h_Uncontained_nc = rt.TH1F("h_Uncontained_nc", "Reco Result of Uncontained Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Range_nc = rt.TH1F("h_Range_nc", "Reco Result of Ranged out Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Flight_decay_nc = rt.TH1F("h_Flight_decay_nc", "Reco Result of Decayed in Flight Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)
h_Sec_int_nc = rt.TH1F("h_Sec_int_nc", "Reco Result of Secondary Interacted Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 5,0,5)

# Distance to Interaction of pions by Reco Status (cc)
h_Dist_reco_nc = rt.TH1F("h_Dist_reco_nc", "Distance to Interaction of Reconstructed Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_misID_nc = rt.TH1F("h_Dist_misID_nc", "Distance to Interaction of mis-ID Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_sec_nc = rt.TH1F("h_Dist_sec_nc", "Distance to Interaction of ID'd as Sec. Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_prim_nc = rt.TH1F("h_Dist_prim_nc", "Distance to Interaction of ID'd as Prim. Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_missing_nc = rt.TH1F("h_Dist_missing_nc", "Distance to Interaction of Missing Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_stack_nc = rt.THStack("h_Dist_stack_nc", "Distance to Interaction of Pions by Reco Type (nc)")
h_Dist_eff_nc = rt.TH1F("h_Dist_eff_nc","Efficiency of Pion Reco by Distance to Interaction (nc);Distance to Interaction (cm);Efficiency(%)",100,0,50)
h_Dist_true_nc = rt.TH1F("h_Dist_true_nc", "Distance to Interaction of True Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_KE_nc = rt.TH2F("h_Dist_KE_nc", "Distance to Interaction vs True Kinetic Energy of Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)

# Distance to Interaction of Pions by History and Reco status (cc)
h_Dist_reco_hist_nc = rt.THStack("h_Dist_reco_hist_nc", "Distance to Interaction of Reconstructed Pions by History (nc)")
h_Dist_misID_hist_nc = rt.THStack("h_Dist_misID_hist_nc", "Distance to Interaction of mis-ID Pions by History (nc)")
h_Dist_sec_hist_nc = rt.THStack("h_Dist_sec_hist_nc", "Distance to Interaction of ID'd as Sec. Pions by History (nc)")
h_Dist_prim_hist_nc = rt.THStack("h_Dist_prim_hist_nc", "Distance to Interaction of ID'd as Prim. Pions by History (nc)")
h_Dist_missing_hist_nc = rt.THStack("h_Dist_missing_hist_nc", "Distance to Interaction of Missing Pions by History (nc)")
Dist_hist_reco_list_nc = [h_Dist_reco_hist_nc,h_Dist_misID_hist_nc,h_Dist_sec_hist_nc,h_Dist_prim_hist_nc,h_Dist_missing_hist_nc]

# KE of easily detectable pions (Colinear angle > 30˚, Range out)
h_range_detect_true_nc = rt.TH1F("h_range_detect_true_nc","True Kinetic Energy of Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_true_cc = rt.TH1F("h_range_detect_true_cc","True Kinetic Energy of Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_reco_nc = rt.TH1F("h_range_detect_reco_nc","True Kinetic Energy of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_reco_cc = rt.TH1F("h_range_detect_reco_cc","True Kinetic Energy of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_eff_cc = rt.TH1F("h_range_detect_eff_cc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by True KE (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_eff_nc = rt.TH1F("h_range_detect_eff_nc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by True KE (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)

detect_count = [0,0]
temp_hists = []

# Temporary histograms for the DOI of pions by history and reco status
for i in range(0,40):
    h_temp = rt.TH1F("","",100,0,100)
    temp_hists.append(h_temp)

# Loop over all entries
for entry in t1:
    p_TID = entry.p_TID
    p_Tracklength = entry.p_Tracklength
    p_TrueE = entry.p_TrueE
    p_Contained = entry.p_Contained
    p_EventID = entry.p_EventID
    p_EventWeight = entry.p_EventWeight
    p_CCNC = entry.p_CCNC
    p_colinearTID = entry.p_colinearTID
    p_colinearAng = entry.p_colinearAng
    p_History = entry.p_History
    p_Dist_to_interaction = entry.p_Dist_to_interaction
    p_Interactions = entry.p_Interactions
    p_Reco_status = entry.p_Reco_status

    nSecondaries = entry.nSecondaries
    sec_TruePDG = entry.sec_TruePDG
    sec_TrueTID = entry.sec_TrueTID
    sec_TrueSource = entry.sec_TrueSource

    KE = p_TrueE[3] - np.sqrt(p_TrueE[3]**2 - (p_TrueE[0]**2 + p_TrueE[1]**2 + p_TrueE[2]**2))

    # CC pions
    if (p_CCNC == 0):
        # Fill relevant truth information
        h_Truecolinear_cc.Fill(p_colinearAng, p_EventWeight)
        h_Dist_true_cc.Fill(p_Dist_to_interaction, p_EventWeight)
        h_Dist_KE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        
        # Easily detectable pions
        if (np.abs(p_colinearAng) < np.sqrt(3)/2):
            if (p_History == 1 and KE > 30):
                h_range_detect_true_cc.Fill(KE, p_EventWeight)
                detect_count[0] += 1
        # Good Recos
        if (p_Reco_status == 0):
            h_Reco_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_cc.Fill(p_History, p_EventWeight)
            h_Dist_reco_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status
            # DOI of good recos by history
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
            
            # Good recos of easily detectable pions
            if (np.abs(p_colinearAng) < np.sqrt(3)/2):
                if (p_History == 1 and KE > 30):
                    h_range_detect_reco_cc.Fill(KE, p_EventWeight)
                    detect_count[1] += 1
        # Mis-IDs
        elif (p_Reco_status == 1):
            h_misID_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_cc.Fill(p_History, p_EventWeight)
            h_Dist_misID_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of mis-IDs by history
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # ID'd as Secondary
        elif (p_Reco_status == 2):
            h_Sec_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_cc.Fill(p_History, p_EventWeight)
            h_Dist_sec_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Sec. IDs by history
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # ID'd as Primary
        elif (p_Reco_status == 3):
            h_Prim_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_prim_cc.Fill(p_History, p_EventWeight)
            h_Dist_prim_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Prim. IDs by history
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # Missing 
        elif (p_Reco_status == 4):
            h_Missing_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_cc.Fill(p_History, p_EventWeight)
            h_Dist_missing_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Missings by history
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)

        # Reco status by history
        if (p_History == 0):
            h_Uncontained_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 1):
            h_Range_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 2):
            h_Flight_decay_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 3):
            h_Sec_int_cc.Fill(p_Reco_status, p_EventWeight)

    # NC pions
    if (p_CCNC == 1):
        # Fill relevant truth information
        h_Truecolinear_nc.Fill(p_colinearAng, p_EventWeight)
        h_Dist_true_nc.Fill(p_Dist_to_interaction, p_EventWeight)
        h_Dist_KE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)

        # Easily reco'd pions
        if (np.abs(p_colinearAng) < np.sqrt(3)/2):
            if (p_History == 1 and KE > 30):
                h_range_detect_true_cc.Fill(KE, p_EventWeight)

        # Good Recos
        if (p_Reco_status == 0):
            h_Reco_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_nc.Fill(p_History, p_EventWeight)
            h_Dist_reco_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of good recos by History
            cat = 4*p_Reco_status + 20
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)

            # good recos of easily reco'd pions
            if (np.abs(p_colinearAng) < np.sqrt(3)/2):
                if (p_History == 1 and KE > 30):
                    h_range_detect_reco_cc.Fill(KE, p_EventWeight)
        # mis-IDs
        elif (p_Reco_status == 1):
            h_misID_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_nc.Fill(p_History, p_EventWeight)
            h_Dist_misID_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of mis-IDs by History
            cat = 4*p_Reco_status + 20
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # ID'd as Secondary
        elif (p_Reco_status == 2):
            h_Sec_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_nc.Fill(p_History, p_EventWeight)
            h_Dist_sec_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Sec. IDs by History
            cat = 4*p_Reco_status + 20
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # ID'd as Primary
        elif (p_Reco_status == 3):
            h_Prim_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_prim_nc.Fill(p_History, p_EventWeight)
            h_Dist_prim_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Prim. IDs by History
            cat = 4*p_Reco_status + 20
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        # Missing
        elif (p_Reco_status == 4):
            h_Missing_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_nc.Fill(p_History, p_EventWeight)
            h_Dist_missing_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            # DOI of Missings by History
            cat = 4*p_Reco_status + 20
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)

        # Reco status by History
        if (p_History == 0):
            h_Uncontained_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 1):
            h_Range_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 2):
            h_Flight_decay_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 3):
            h_Sec_int_nc.Fill(p_Reco_status, p_EventWeight)

# Scale histograms to target POT
h_Reco_colinear_cc.Scale(6.67e20/t_pot)
h_misID_colinear_cc.Scale(6.67e20/t_pot)
h_Sec_colinear_cc.Scale(6.67e20/t_pot)
h_Prim_colinear_cc.Scale(6.67e20/t_pot)
h_Missing_colinear_cc.Scale(6.67e20/t_pot)
h_Truecolinear_cc.Scale(6.67e20/t_pot)
h_Hist_reco_cc.Scale(6.67e20/t_pot)
h_Hist_misID_cc.Scale(6.67e20/t_pot)
h_Hist_sec_cc.Scale(6.67e20/t_pot)
h_Hist_prim_cc.Scale(6.67e20/t_pot)
h_Hist_missing_cc.Scale(6.67e20/t_pot)
h_Uncontained_cc.Scale(6.67e20/t_pot)
h_Range_cc.Scale(6.67e20/t_pot)
h_Flight_decay_cc.Scale(6.67e20/t_pot)
h_Sec_int_cc.Scale(6.67e20/t_pot)

h_Dist_reco_cc.Scale(6.67e20/t_pot)
h_Dist_misID_cc.Scale(6.67e20/t_pot) 
h_Dist_sec_cc.Scale(6.67e20/t_pot) 
h_Dist_prim_cc.Scale(6.67e20/t_pot) 
h_Dist_missing_cc.Scale(6.67e20/t_pot)
h_Dist_true_cc.Scale(6.67e20/t_pot)
h_Effcolinear_cc.Divide(h_Reco_colinear_cc, h_Truecolinear_cc)
h_Dist_eff_cc.Divide(h_Dist_reco_cc, h_Dist_true_cc)

h_Reco_colinear_nc.Scale(6.67e20/t_pot)
h_misID_colinear_nc.Scale(6.67e20/t_pot)
h_Sec_colinear_nc.Scale(6.67e20/t_pot)
h_Missing_colinear_nc.Scale(6.67e20/t_pot)
h_Truecolinear_nc.Scale(6.67e20/t_pot)
h_Hist_reco_nc.Scale(6.67e20/t_pot)
h_Hist_misID_nc.Scale(6.67e20/t_pot)
h_Hist_sec_nc.Scale(6.67e20/t_pot)
h_Hist_missing_nc.Scale(6.67e20/t_pot)
h_Uncontained_nc.Scale(6.67e20/t_pot)
h_Range_nc.Scale(6.67e20/t_pot)
h_Flight_decay_nc.Scale(6.67e20/t_pot)
h_Sec_int_nc.Scale(6.67e20/t_pot)

h_Dist_reco_nc.Scale(6.67e20/t_pot)
h_Dist_misID_nc.Scale(6.67e20/t_pot) 
h_Dist_sec_nc.Scale(6.67e20/t_pot) 
h_Dist_prim_nc.Scale(6.67e20/t_pot) 
h_Dist_missing_nc.Scale(6.67e20/t_pot)
h_Dist_true_nc.Scale(6.67e20/t_pot)
h_Effcolinear_nc.Divide(h_Reco_colinear_nc, h_Truecolinear_nc)
h_Dist_eff_nc.Divide(h_Dist_reco_nc, h_Dist_true_nc)

for i in range(0,40):
    temp_hists[i].Scale(6.67e20/t_pot)
    print(temp_hists[i].GetEntries())

h_range_detect_reco_cc.Scale(6.67e20/t_pot)
h_range_detect_reco_nc.Scale(6.67e20/t_pot)
h_range_detect_true_cc.Scale(6.67e20/t_pot)
h_range_detect_true_nc.Scale(6.67e20/t_pot)
h_range_detect_eff_cc.Divide(h_range_detect_reco_cc,h_range_detect_true_cc)
h_range_detect_eff_nc.Divide(h_range_detect_reco_nc,h_range_detect_true_nc)

# Stack histogram function
History_tags = ["Uncontained", "Range Out", "Sec. Interaction", "Decay in Flight"]
Reco_tags = ["Good Reco", "Mis-ID", "ID'd as Sec.", "ID'd as Prim.", "Missing"]
def write_stack(h_list, h_stack, name, l_pos, l_tags):
    canvas = rt.TCanvas(str(name))
    colors = [rt.kGreen-3,rt.kBlue-7,rt.kYellow-7,rt.kViolet-4,rt.kRed-7]
    for i in range(0,len(h_list)-1):
        h_temp = h_list[0].Clone()
        h_temp.SetFillColor(colors[i])
        h_temp.Draw("HIST")
        h_stack.Add(h_temp)
    h_temp = h_list[-1].Clone()
    h_temp.SetFillColor(colors[4])
    legend = rt.TLegend(l_pos[0], l_pos[1], l_pos[2], l_pos[3])
    for i in range(0,len(h_list)-1):
        legend.AddEntry(h_list[i],l_tags[1],"f")
    legend.AddEntry(h_temp,l_tags[-1],"f")
    h_temp.GetListOfFunctions().Add(legend)
    h_temp.Draw()
    legend.Draw()
    h_stack.Add(h_temp)
    h_stack.Draw("HIST")
    return canvas

f_out.cd()
# Create and write stacked histograms
# Collinearity, DOI stacks
c_Stack_colinear_cc = write_stack([h_Reco_colinear_cc,h_misID_colinear_cc,
                                  h_Sec_colinear_cc,h_Prim_colinear_cc,h_Missing_colinear_cc],
                                  h_Stack_colinear_cc,"c_Stack_colinear_cc",
                                  [0.1, 0.5, 0.3, 0.9], Reco_tags)
h_Stack_colinear_cc.GetXaxis().SetTitle("Most Colinear Primary Angle (cos theta)")
h_Stack_colinear_cc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Stack_colinear_cc.Update()
c_Stack_colinear_cc.Write()

c_Dist_stack_cc = write_stack([h_Dist_reco_cc,h_Dist_misID_cc,h_Dist_sec_cc,h_Dist_prim_cc,
                              h_Dist_missing_cc],h_Dist_stack_cc,"c_Dist_stack_cc",
                              [0.7, 0.5, 0.9, 0.9], Reco_tags)
h_Dist_stack_cc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_cc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_cc.Update()
c_Dist_stack_cc.Write()

c_Stack_colinear_nc = write_stack([h_Reco_colinear_nc,h_misID_colinear_nc,
                                  h_Sec_colinear_nc,h_Prim_colinear_nc,h_Missing_colinear_nc],
                                  h_Stack_colinear_nc,"c_Stack_colinear_nc",
                                  [0.1, 0.5, 0.3, 0.9], Reco_tags)
h_Stack_colinear_nc.GetXaxis().SetTitle("Most Colinear Primary Angle (cos theta)")
h_Stack_colinear_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Stack_colinear_nc.Update()
c_Stack_colinear_nc.Write()

c_Dist_stack_nc = write_stack([h_Dist_reco_nc,h_Dist_misID_nc,h_Dist_sec_nc,h_Dist_prim_nc,
                              h_Dist_missing_nc],h_Dist_stack_nc,"c_Dist_stack_nc",
                              [0.7, 0.5, 0.9, 0.9], Reco_tags)
h_Dist_stack_nc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_nc.Update()
c_Dist_stack_nc.Write()

# DOI by history and reco status stacks
nc_stack_names = ["c_Dist_reco_hist_nc","c_Dist_misID_hist_nc","c_Dist_sec_hist_nc","c_Dist_prim_hist_nc","c_Dist_missing_hist_nc"]
for i in range(len(nc_stack_names),10):
    canvas = write_stack([temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3]],
                Dist_hist_reco_list_nc[i-5],nc_stack_names[i-5],[0.7, 0.55, 0.9, 0.9], History_tags)
    Dist_hist_reco_list_nc[i-5].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_nc[i-5].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

cc_stack_names = ["c_Dist_reco_hist_cc","c_Dist_misID_hist_cc","c_Dist_sec_hist_cc","c_Dist_prim_hist_cc","c_Dist_missing_hist_cc"]
for i in range(0,len(cc_stack_names)):
    canvas = write_stack([temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3]],
                Dist_hist_reco_list_cc[i],cc_stack_names[i],[0.7, 0.55, 0.9, 0.9], History_tags)
    Dist_hist_reco_list_cc[i].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_cc[i].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

# Draw histograms and fill bin tags
h_Reco_colinear_cc.Draw()
h_misID_colinear_cc.Draw()
h_Sec_colinear_cc.Draw()
h_Missing_colinear_cc.Draw()
h_Truecolinear_cc.Draw()
h_Effcolinear_cc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_reco_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_misID_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_sec_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_prim_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_missing_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_reco_cc.Draw()
h_Hist_misID_cc.Draw()
h_Hist_sec_cc.Draw()
h_Hist_prim_cc.Draw()
h_Hist_missing_cc.Draw()

for i in range(0, len(Reco_tags)):
    h_Uncontained_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Range_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Flight_decay_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Sec_int_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Uncontained_cc.Draw()
h_Range_cc.Draw()
h_Flight_decay_cc.Draw()
h_Sec_int_cc.Draw()

h_Reco_colinear_nc.Draw()
h_misID_colinear_nc.Draw()
h_Sec_colinear_nc.Draw()
h_Missing_colinear_nc.Draw()
h_Truecolinear_nc.Draw()
h_Effcolinear_nc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_reco_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_misID_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_sec_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_prim_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
    h_Hist_missing_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_reco_nc.Draw()
h_Hist_misID_nc.Draw()
h_Hist_sec_nc.Draw()
h_Hist_prim_nc.Draw()
h_Hist_missing_nc.Draw()

for i in range(0, len(Reco_tags)):
    h_Uncontained_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Range_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Flight_decay_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
    h_Sec_int_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Uncontained_nc.Draw()
h_Range_nc.Draw()
h_Flight_decay_nc.Draw()
h_Sec_int_nc.Draw()

# Single histogram function
def write_canvas(h1, name):
    c_name = "c" + name[1:]
    canvas = rt.TCanvas(str(c_name))
    # h_temp.SetFillColor(rt.kViolet-4)
    h1.Draw("HIST")
    canvas.Update()
    canvas.Draw("HIST")
    canvas.Write()
# write_canvas(h_Reco_colinear_cc, "h_Reco_colinear_cc")
# write_canvas(h_misID_colinear_cc, "h_misID_colinear_cc")
# write_canvas(h_Sec_colinear_cc, "h_Sec_colinear_cc")
write_canvas(h_Prim_colinear_cc, "h_Prim_colinear_cc")
# write_canvas(h_Missing_colinear_cc, "h_Missing_colinear_cc")
# write_canvas(h_Truecolinear_cc, "h_Truecolinear_cc")
write_canvas(h_Effcolinear_cc, "h_Effcolinear_cc")

write_canvas(h_Hist_reco_cc, "h_Hist_reco_cc")
write_canvas(h_Hist_misID_cc, "h_Hist_misID_cc")
write_canvas(h_Hist_sec_cc, "h_Hist_sec_cc")
write_canvas(h_Hist_prim_cc, "h_Hist_prim_cc")
write_canvas(h_Hist_missing_cc, "h_Hist_missing_cc")

write_canvas(h_Uncontained_cc, "h_Uncontained_cc")
write_canvas(h_Range_cc, "h_Range_cc")
write_canvas(h_Flight_decay_cc, "h_Flight_decay_cc")
write_canvas(h_Sec_int_cc, "h_Sec_int_cc")


# write_canvas(h_Reco_colinear_nc, "h_Reco_colinear_nc")
# write_canvas(h_misID_colinear_nc, "h_misID_colinear_nc")
# write_canvas(h_Sec_colinear_nc, "h_Sec_colinear_nc")
write_canvas(h_Prim_colinear_nc, "h_prim_colinear_nc")
# write_canvas(h_Missing_colinear_nc, "h_Missing_colinear_nc")
# write_canvas(h_Truecolinear_nc, "h_Truecolinear_nc")
write_canvas(h_Effcolinear_nc, "h_Effcolinear_nc")

write_canvas(h_Hist_reco_nc, "h_Hist_reco_nc")
write_canvas(h_Hist_misID_nc, "h_Hist_misID_nc")
write_canvas(h_Hist_sec_nc, "h_Hist_sec_nc")
write_canvas(h_Hist_prim_nc, "h_Hist_prim_nc")
write_canvas(h_Hist_missing_nc, "h_Hist_missing_nc")

write_canvas(h_Uncontained_nc, "h_Uncontained_nc")
write_canvas(h_Range_nc, "h_Range_nc")
write_canvas(h_Flight_decay_nc, "h_Flight_decay_nc")
write_canvas(h_Sec_int_nc, "h_Sec_int_nc")

write_canvas(h_Dist_eff_cc, "h_Dist_eff_cc")
write_canvas(h_Dist_eff_nc, "h_Dist_eff_nc")

write_canvas(h_range_detect_eff_cc, "h_range_detect_eff_cc")
write_canvas(h_range_detect_eff_nc, "h_range_detect_eff_nc")

# 2D histogram of KE vs. DOI
c_Dist_KE_cc = rt.TCanvas("c_Dist_KE_cc")
h_Dist_KE_cc.Draw("Colz")
c_Dist_KE_cc.Update()
c_Dist_KE_cc.Draw()
c_Dist_KE_cc.Write()
c_Dist_KE_nc = rt.TCanvas("c_Dist_KE_nc")
h_Dist_KE_nc.Draw("Colz")
c_Dist_KE_nc.Update()
c_Dist_KE_nc.Draw()
c_Dist_KE_nc.Write()
f_out.Close()


# print(detect_count[0], detect_count[1],detect_count[1]/detect_count[0])