# dist_eff.py
# Author: Andy Zhang
# Purpose: Find Reconstruction Efficiencies of True Primary Pions ordered by 
# distance to first interaction. Counterpart of hist_eff.py
# Data Created: July 31st, 2024
# Last Modified: July 31st, 2024

import ROOT as rt
import numpy as np
from ROOT import gROOT, AddressOf
f1 = rt.TFile("pion_stats.root")
f_out = rt.TFile("dist_eff.root","RECREATE")
t1 = f1.Get("pion_tree")
p1 = f1.Get("potTree")

# Count total good POT for scaling
t_pot = 0.
for entry in p1:
    t_pot += entry.totGoodPOT

# Initialize all histograms

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
h_Stack_DKE_cc = rt.THStack("h_Stack_DKE_cc", "Distance to Interaction vs True Kinetic Energy of Pions (cc)")
h_Uncontained_DKE_cc = rt.TH2F("h_Uncontained_DKE_cc", "Distance to Interaction vs True Kinetic Energy of Uncontained Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Range_DKE_cc = rt.TH2F("h_Range_DKE_cc", "Distance to Interaction vs True Kinetic Energy of Ranged Out Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Sec_DKE_cc = rt.TH2F("h_Sec_DKE_cc", "Distance to Interaction vs True Kinetic Energy of Pions with Sec. Int.;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Decay_DKE_cc = rt.TH2F("h_Decay_DKE_cc", "Distance to Interaction vs True Kinetic Energy of Decayed in Flight Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)

h_Dist_ang_true_cc = rt.TH2F("h_Dist_ang_true_cc", "Distance to Interaction vs. Nearest Collinear Primary of True Pions;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)
h_Dist_ang_reco_cc = rt.TH2F("h_Dist_ang_reco_cc", "Distance to Interaction vs. Nearest Collinear Primary of Reco Pions;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)
h_Dist_ang_eff_cc = rt.TH2F("h_Dist_ang_eff_cc", "Pion Reco Efficiency by Distance to Interaction vs. Nearest Collinear Primary;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)

# Distance to Interaction of Pions by History and Reco status (cc)
h_Dist_reco_hist_cc = rt.THStack("h_Dist_reco_hist_cc", "Distacce to Interaction of Reconstructed Pions by History (cc)")
h_Dist_misID_hist_cc = rt.THStack("h_Dist_misID_hist_cc", "Distacce to Interaction of mis-ID Pions by History (cc)")
h_Dist_sec_hist_cc = rt.THStack("h_Dist_sec_hist_cc", "Distacce to Interaction of ID'd as Sec. Pions by History (cc)")
h_Dist_prim_hist_cc = rt.THStack("h_Dist_prim_hist_cc", "Distance to Interaction of ID'd as Prim. Pions by History (cc)")
h_Dist_missing_hist_cc = rt.THStack("h_Dist_missing_hist_cc", "Distacce to Interaction of Missing Pions by History (cc)")
Dist_hist_reco_list_cc = [h_Dist_reco_hist_cc,h_Dist_misID_hist_cc,h_Dist_sec_hist_cc,h_Dist_prim_hist_cc,h_Dist_missing_hist_cc]

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
h_Stack_DKE_nc = rt.THStack("h_Stack_DKE_nc", "Distance to Interaction vs True Kinetic Energy of Pions (nc)")
h_Uncontained_DKE_nc = rt.TH2F("h_Uncontained_DKE_nc", "Distance to Interaction vs True Kinetic Energy of Uncontained Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Range_DKE_nc = rt.TH2F("h_Range_DKE_nc", "Distance to Interaction vs True Kinetic Energy of Ranged Out Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Sec_DKE_nc = rt.TH2F("h_Sec_DKE_nc", "Distance to Interaction vs True Kinetic Energy of Pions with Sec. Int.;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)
h_Decay_DKE_nc = rt.TH2F("h_Decay_DKE_nc", "Distance to Interaction vs True Kinetic Energy of Decayed in Flight Pions;True Kinetic Energy (MeV);Distance to Interaction (cm)",100,0,600,150,0,150)

h_Dist_ang_true_nc = rt.TH2F("h_Dist_ang_true_nc", "Distance to Interaction vs. Nearest Collinear Primary of True Pions;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)
h_Dist_ang_reco_nc = rt.TH2F("h_Dist_ang_reco_nc", "Distance to Interaction vs. Nearest Collinear Primary of Reco Pions;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)
h_Dist_ang_eff_nc = rt.TH2F("h_Dist_ang_eff_nc", "Pion Reco Efficiency by Distance to Interaction vs. Nearest Collinear Primary;Nearest Collinear Angle (cos theta);Distance to Interaction (cm)", 40,-1,1,150,0,150)

# Distance to Interaction of Pions by History and Reco status (nc)
h_Dist_reco_hist_nc = rt.THStack("h_Dist_reco_hist_nc", "Distance to Interaction of Reconstructed Pions by History (nc)")
h_Dist_misID_hist_nc = rt.THStack("h_Dist_misID_hist_nc", "Distance to Interaction of mis-ID Pions by History (nc)")
h_Dist_sec_hist_nc = rt.THStack("h_Dist_sec_hist_nc", "Distance to Interaction of ID'd as Sec. Pions by History (nc)")
h_Dist_prim_hist_nc = rt.THStack("h_Dist_prim_hist_nc", "Distance to Interaction of ID'd as Prim. Pions by History (nc)")
h_Dist_missing_hist_nc = rt.THStack("h_Dist_missing_hist_nc", "Distance to Interaction of Missing Pions by History (nc)")
Dist_hist_reco_list_nc = [h_Dist_reco_hist_nc,h_Dist_misID_hist_nc,h_Dist_sec_hist_nc,h_Dist_prim_hist_nc,h_Dist_missing_hist_nc]

# KE of easily detectable pions (Colinear angle > 30˚, Range out)
h_range_detect_true_nc = rt.TH1F("h_range_detect_true_nc","Distance to Interaction of Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",60,0,120)
h_range_detect_true_cc = rt.TH1F("h_range_detect_true_cc","Distance to Interaction of Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",60,0,120)
h_range_detect_reco_nc = rt.TH1F("h_range_detect_reco_nc","Distance to Interaction of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",60,0,120)
h_range_detect_reco_cc = rt.TH1F("h_range_detect_reco_cc","Distance to Interaction of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",60,0,120)
h_range_detect_eff_cc = rt.TH1F("h_range_detect_eff_cc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by DOI (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT",60,0,120)
h_range_detect_eff_nc = rt.TH1F("h_range_detect_eff_nc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by DOI (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT",60,0,120)


detect_count = [0,0]
temp_hists = []
prim_counter = 0
# Temporary histograms for the DOI of pions by history and reco status
for i in range(0,40):
    h_temp = rt.TH1F("","",100,0,100)
    temp_hists.append(h_temp)

missed_long_pions = []
missed_long_pions_counter = 0
# Loop over all entries
for entry in t1:
    p_TID = entry.p_TID
    p_Tracklength = entry.p_Tracklength
    p_TrueE = entry.p_TrueE
    p_Contained = entry.p_Contained

    p_FileID = entry.p_FileID
    p_EventID = entry.p_EventID
    p_Run = entry.p_Run
    p_Subrun = entry.p_Subrun
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

    if (p_Dist_to_interaction > 15 and p_Reco_status == 4):
        pion_event = [p_FileID,p_Run,p_Subrun,p_EventID]
        missed_long_pions.append(pion_event)
        missed_long_pions_counter += 1

    KE = p_TrueE[3] - np.sqrt(p_TrueE[3]**2 - (p_TrueE[0]**2 + p_TrueE[1]**2 + p_TrueE[2]**2))

    # CC pions
    if (p_CCNC == 0):
        # Fill relevant truth information
        h_Dist_true_cc.Fill(p_Dist_to_interaction, p_EventWeight)
        h_Dist_KE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        h_Dist_ang_true_cc.Fill(p_colinearAng, p_Dist_to_interaction, p_EventWeight)
        if (np.abs(p_colinearAng) < np.sqrt(3)/2 and p_History == 1):
            h_range_detect_true_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            detect_count[0] += 1
        # Good Recos
        if (p_Reco_status == 0):
            h_Dist_reco_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            h_Dist_ang_reco_cc.Fill(p_colinearAng, p_Dist_to_interaction, p_EventWeight)
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
            if (np.abs(p_colinearAng) < np.sqrt(3)/2 and p_History == 1):
                h_range_detect_reco_cc.Fill(p_Dist_to_interaction, p_EventWeight)
                detect_count[1] += 1
            
        # Mis-IDs
        elif (p_Reco_status == 1):
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

        if (p_History == 0):
            h_Uncontained_DKE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 1):
            h_Range_DKE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 2):
            h_Sec_DKE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 3):
            h_Decay_DKE_cc.Fill(KE, p_Dist_to_interaction, p_EventWeight)

    # NC pions
    if (p_CCNC == 1):
        # Fill relevant truth information
        h_Dist_true_nc.Fill(p_Dist_to_interaction, p_EventWeight)
        h_Dist_KE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        h_Dist_ang_true_nc.Fill(p_colinearAng, p_Dist_to_interaction, p_EventWeight)

        if (np.abs(p_colinearAng) < np.sqrt(3)/2 and p_History == 1):
            h_range_detect_true_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            detect_count[0] += 1
        # Good Recos
        if (p_Reco_status == 0):
            h_Dist_reco_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            h_Dist_ang_reco_nc.Fill(p_colinearAng, p_Dist_to_interaction, p_EventWeight)
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
        if (np.abs(p_colinearAng) < np.sqrt(3)/2 and p_History == 1):
            h_range_detect_reco_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            detect_count[1] += 1

        # mis-IDs
        elif (p_Reco_status == 1):
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
            prim_counter += 1
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

        if (p_History == 0):
            h_Uncontained_DKE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 1):
            h_Range_DKE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 2):
            h_Sec_DKE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)
        elif (p_History == 3):
            h_Decay_DKE_nc.Fill(KE, p_Dist_to_interaction, p_EventWeight)

# Scale histograms to target POT
h_Dist_reco_cc.Scale(6.67e20/t_pot)
h_Dist_misID_cc.Scale(6.67e20/t_pot) 
h_Dist_sec_cc.Scale(6.67e20/t_pot) 
h_Dist_prim_cc.Scale(6.67e20/t_pot) 
h_Dist_missing_cc.Scale(6.67e20/t_pot)
h_Dist_true_cc.Scale(6.67e20/t_pot)
h_Dist_ang_true_cc.Scale(6.67e20/t_pot)
h_Dist_ang_reco_cc.Scale(6.67e20/t_pot)
h_Dist_eff_cc.Divide(h_Dist_reco_cc, h_Dist_true_cc)
h_Dist_ang_eff_cc.Divide(h_Dist_ang_reco_cc, h_Dist_ang_true_cc)

h_Dist_reco_nc.Scale(6.67e20/t_pot)
h_Dist_misID_nc.Scale(6.67e20/t_pot) 
h_Dist_sec_nc.Scale(6.67e20/t_pot)
h_Dist_prim_nc.Scale(6.67e20/t_pot)
h_Dist_missing_nc.Scale(6.67e20/t_pot)
h_Dist_true_nc.Scale(6.67e20/t_pot)
h_Dist_ang_true_nc.Scale(6.67e20/t_pot)
h_Dist_ang_reco_nc.Scale(6.67e20/t_pot)
h_Dist_eff_nc.Divide(h_Dist_reco_nc, h_Dist_true_nc)
h_Dist_ang_eff_nc.Divide(h_Dist_ang_reco_nc, h_Dist_ang_true_nc)

h_range_detect_reco_cc.Scale(6.67e20/t_pot)
h_range_detect_reco_nc.Scale(6.67e20/t_pot)
h_range_detect_true_cc.Scale(6.67e20/t_pot)
h_range_detect_true_nc.Scale(6.67e20/t_pot)
h_range_detect_eff_cc.Divide(h_range_detect_reco_cc,h_range_detect_true_cc)
h_range_detect_eff_nc.Divide(h_range_detect_reco_nc,h_range_detect_true_nc)


for i in range(0,40):
    temp_hists[i].Scale(6.67e20/t_pot)
    # print(temp_hists[i].GetEntries())

# Stack histogram function
History_tags = ["Range Out", "Sec. Interaction", "Decay in Flight", "Uncontained"]
Reco_tags = ["Good Reco", "Mis-ID", "ID'd as Sec.", "ID'd as Prim.", "Missing"]
def write_stack(h_list, h_stack, name, l_pos, l_tags, draw_option):
    canvas = rt.TCanvas(str(name))
    colors = [rt.kGreen-3,rt.kBlue-7,rt.kYellow-7,rt.kViolet-4,rt.kRed-7]
    h_temp = [rt.TH1F(),rt.TH1F(),rt.TH1F(),rt.TH1F(),rt.TH1F()]
    for i in range(0,len(h_list)-1):
        h_temp[i] = h_list[i].Clone()
        h_temp[i].SetFillColor(colors[i])
        h_temp[i].Draw(str(draw_option))
        h_stack.Add(h_temp[i])
    h_temp[-1] = h_list[-1].Clone()
    h_temp[-1].SetFillColor(colors[4])
    legend = rt.TLegend(l_pos[0], l_pos[1], l_pos[2], l_pos[3])
    for i in range(0,len(h_list)-1):
        legend.AddEntry(h_list[i],l_tags[1],"f")
    legend.AddEntry(h_temp[-1],l_tags[-1],"f")
    h_temp[-1].GetListOfFunctions().Add(legend)
    h_temp[-1].Draw(str(draw_option))
    legend.Draw(str(draw_option))
    h_stack.Add(h_temp[-1])
    h_stack.Draw(str(draw_option))
    return canvas

f_out.cd()
# Create and write stacked histograms
# DOI stacks
c_Dist_stack_cc = write_stack([h_Dist_reco_cc,h_Dist_misID_cc,h_Dist_sec_cc,h_Dist_prim_cc,
                              h_Dist_missing_cc],h_Dist_stack_cc,"c_Dist_stack_cc",
                              [0.7, 0.5, 0.9, 0.9], Reco_tags, "HIST")
h_Dist_stack_cc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_cc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_cc.Update()
c_Dist_stack_cc.Write()

c_Dist_stack_nc = write_stack([h_Dist_reco_nc,h_Dist_misID_nc,h_Dist_sec_nc,h_Dist_prim_nc,
                              h_Dist_missing_nc],h_Dist_stack_nc,"c_Dist_stack_nc",
                              [0.7, 0.5, 0.9, 0.9], Reco_tags, "HIST")
h_Dist_stack_nc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_nc.Update()
c_Dist_stack_nc.Write()

# DOI by history and reco status stacks
nc_stack_names = ["c_Dist_reco_hist_nc","c_Dist_misID_hist_nc","c_Dist_sec_hist_nc","c_Dist_prim_hist_nc","c_Dist_missing_hist_nc"]
for i in range(len(nc_stack_names),10):
    canvas = write_stack([temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3]],
                Dist_hist_reco_list_nc[i-5],nc_stack_names[i-5],[0.7, 0.55, 0.9, 0.9], History_tags, "HIST")
    Dist_hist_reco_list_nc[i-5].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_nc[i-5].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

cc_stack_names = ["c_Dist_reco_hist_cc","c_Dist_misID_hist_cc","c_Dist_sec_hist_cc","c_Dist_prim_hist_cc","c_Dist_missing_hist_cc"]
for i in range(0,len(cc_stack_names)):
    canvas = write_stack([temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3]],
                Dist_hist_reco_list_cc[i],cc_stack_names[i],[0.7, 0.55, 0.9, 0.9], History_tags, "HIST")
    Dist_hist_reco_list_cc[i].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_cc[i].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

c_Stack_DKE_cc = write_stack([h_Range_DKE_cc,h_Sec_DKE_cc,h_Decay_DKE_cc,h_Uncontained_DKE_cc],
                             h_Stack_DKE_cc,"c_Stack_DKE_cc",
                             [0.7, 0.5, 0.9, 0.9], History_tags, "LEGO1")
h_Stack_DKE_cc.GetXaxis().SetTitle("True Pion Kinetic Energy (MeV)")
h_Stack_DKE_cc.GetYaxis().SetTitle("Distance to Interaction (cm)")
c_Stack_DKE_cc.Update()
c_Stack_DKE_cc.Write()

c_Stack_DKE_nc = write_stack([h_Range_DKE_nc,h_Sec_DKE_nc,h_Decay_DKE_nc,h_Uncontained_DKE_nc],
                             h_Stack_DKE_nc,"c_Stack_DKE_nc",
                             [0.7, 0.5, 0.9, 0.9], History_tags, "LEGO1")
h_Stack_DKE_nc.GetXaxis().SetTitle("True Pion Kinetic Energy (MeV)")
h_Stack_DKE_nc.GetYaxis().SetTitle("Distance to Interaction (cm)")
c_Stack_DKE_nc.Update()
c_Stack_DKE_nc.Write()

# Single histogram function
def write_canvas(h1, name, draw_option):
    c_name = "c" + name[1:]
    canvas = rt.TCanvas(str(c_name))
    # h_temp.SetFillColor(rt.kViolet-4)
    h1.Draw(str(draw_option))
    canvas.Update()
    canvas.Draw(str(draw_option))
    canvas.Write()

write_canvas(h_Dist_eff_cc, "h_Dist_eff_cc", "HIST")
write_canvas(h_Dist_eff_nc, "h_Dist_eff_nc", "HIST")

# 2D histogram of KE vs. DOI
write_canvas(h_Dist_KE_cc, "h_Dist_KE_cc", "LEGO1")
write_canvas(h_Dist_KE_nc, "h_Dist_KE_nc", "LEGO1")
write_canvas(h_Uncontained_DKE_cc, "h_Uncontained_DKE_cc", "LEGO1")
write_canvas(h_Uncontained_DKE_nc, "h_Uncontained_DKE_nc", "LEGO1")
write_canvas(h_Range_DKE_cc, "h_Range_DKE_cc", "LEGO1")
write_canvas(h_Range_DKE_nc, "h_Range_DKE_nc", "LEGO1")
write_canvas(h_Sec_DKE_cc, "h_Sec_DKE_cc", "LEGO1")
write_canvas(h_Sec_DKE_nc, "h_Sec_DKE_nc", "LEGO1")
write_canvas(h_Decay_DKE_cc, "h_Decay_DKE_cc", "LEGO1")
write_canvas(h_Decay_DKE_nc, "h_Decay_DKE_nc", "LEGO1")

write_canvas(h_range_detect_eff_cc, "h_range_detect_eff_cc", "HIST")
write_canvas(h_range_detect_eff_nc, "h_range_detect_eff_nc", "HIST")

write_canvas(h_Dist_ang_eff_cc, "h_Dist_ang_eff_cc", "COLZ")
write_canvas(h_Dist_ang_eff_nc, "h_Dist_ang_eff_nc", "COLZ")
f_out.Close()

print(detect_count[0], detect_count[1],detect_count[1]/detect_count[0])