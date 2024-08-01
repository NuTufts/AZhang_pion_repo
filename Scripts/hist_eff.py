# hist_eff.py
# Author: Andy Zhang
# Purpose: Find Reconstruction Efficiencies of True Primary Pions ordered by pion history,
# reco status type, and most colinear primary. Counterpart of dist_eff.py
# Data Created: June 5th, 2024
# Last Modified: July 31st, 2024

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

# KE of easily detectable pions (Colinear angle > 30˚, Range out)
h_range_detect_true_nc = rt.TH1F("h_range_detect_true_nc","True Kinetic Energy of Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_true_cc = rt.TH1F("h_range_detect_true_cc","True Kinetic Energy of Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_reco_nc = rt.TH1F("h_range_detect_reco_nc","True Kinetic Energy of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_reco_cc = rt.TH1F("h_range_detect_reco_cc","True Kinetic Energy of Reco'd Range Out Pions with Max Colinear Angle > 30˚ (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_eff_cc = rt.TH1F("h_range_detect_eff_cc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by True KE (cc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)
h_range_detect_eff_nc = rt.TH1F("h_range_detect_eff_nc","Efficiency of Range Out Pions with Max Colinear Angle > 30˚ by True KE (nc);Kinetic Energy (MeV);# of Pions per 6.67e20 POT",40,0,1200)

detect_count = [0,0]
prim_counter = 0
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

    if (p_Dist_to_interaction > 15 and p_Reco_status == 4 and p_colinearAng < np.sqrt(3)/2):
        pion_event = [p_FileID,p_Run,p_Subrun,p_EventID]
        missed_long_pions.append(pion_event)
        missed_long_pions_counter += 1

    KE = p_TrueE[3] - np.sqrt(p_TrueE[3]**2 - (p_TrueE[0]**2 + p_TrueE[1]**2 + p_TrueE[2]**2))

    # CC pions
    if (p_CCNC == 0):
        # Fill relevant truth information
        h_Truecolinear_cc.Fill(p_colinearAng, p_EventWeight)

        # Easily detectable pions
        if (np.abs(p_colinearAng) < np.sqrt(3)/2):
            if (p_History == 1 and KE > 30):
                h_range_detect_true_cc.Fill(KE, p_EventWeight)
                detect_count[0] += 1
        # Good Recos
        if (p_Reco_status == 0):
            h_Reco_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_cc.Fill(p_History, p_EventWeight)
            cat = 4*p_Reco_status
            
            # Good recos of easily detectable pions
            if (np.abs(p_colinearAng) < np.sqrt(3)/2):
                if (p_History == 1 and KE > 30):
                    h_range_detect_reco_cc.Fill(KE, p_EventWeight)
                    detect_count[1] += 1
        # Mis-IDs
        elif (p_Reco_status == 1):
            h_misID_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_cc.Fill(p_History, p_EventWeight)
        # ID'd as Secondary
        elif (p_Reco_status == 2):
            h_Sec_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_cc.Fill(p_History, p_EventWeight)
        # ID'd as Primary
        elif (p_Reco_status == 3):
            h_Prim_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_prim_cc.Fill(p_History, p_EventWeight)
        # Missing 
        elif (p_Reco_status == 4):
            h_Missing_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_cc.Fill(p_History, p_EventWeight)

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

        # Easily reco'd pions
        if (np.abs(p_colinearAng) < np.sqrt(3)/2):
            if (p_History == 1 and KE > 30):
                h_range_detect_true_cc.Fill(KE, p_EventWeight)

        # Good Recos
        if (p_Reco_status == 0):
            h_Reco_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_nc.Fill(p_History, p_EventWeight)

            # good recos of easily reco'd pions
            if (np.abs(p_colinearAng) < np.sqrt(3)/2):
                if (p_History == 1 and KE > 30):
                    h_range_detect_reco_cc.Fill(KE, p_EventWeight)
        # mis-IDs
        elif (p_Reco_status == 1):
            h_misID_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_nc.Fill(p_History, p_EventWeight)
        # ID'd as Secondary
        elif (p_Reco_status == 2):
            h_Sec_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_nc.Fill(p_History, p_EventWeight)
        # ID'd as Primary
        elif (p_Reco_status == 3):
            prim_counter += 1
            h_Prim_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_prim_nc.Fill(p_History, p_EventWeight)
        # Missing
        elif (p_Reco_status == 4):
            h_Missing_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_nc.Fill(p_History, p_EventWeight)

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
h_Effcolinear_cc.Divide(h_Reco_colinear_cc, h_Truecolinear_cc)

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
h_Effcolinear_nc.Divide(h_Reco_colinear_nc, h_Truecolinear_nc)

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
    h_temp = [rt.TH1F(),rt.TH1F(),rt.TH1F(),rt.TH1F(),rt.TH1F()]
    for i in range(0,len(h_list)-1):
        h_temp[i] = h_list[i].Clone()
        h_temp[i].SetFillColor(colors[i])
        h_temp[i].Draw("HIST")
        h_stack.Add(h_temp[i])
    h_temp[-1] = h_list[-1].Clone()
    h_temp[-1].SetFillColor(colors[4])
    legend = rt.TLegend(l_pos[0], l_pos[1], l_pos[2], l_pos[3])
    for i in range(0,len(h_list)-1):
        legend.AddEntry(h_list[i],l_tags[1],"f")
    legend.AddEntry(h_temp[-1],l_tags[-1],"f")
    h_temp[-1].GetListOfFunctions().Add(legend)
    h_temp[-1].Draw()
    legend.Draw()
    h_stack.Add(h_temp[-1])
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

c_Stack_colinear_nc = write_stack([h_Reco_colinear_nc,h_misID_colinear_nc,
                                  h_Sec_colinear_nc,h_Prim_colinear_nc,h_Missing_colinear_nc],
                                  h_Stack_colinear_nc,"c_Stack_colinear_nc",
                                  [0.1, 0.5, 0.3, 0.9], Reco_tags)
h_Stack_colinear_nc.GetXaxis().SetTitle("Most Colinear Primary Angle (cos theta)")
h_Stack_colinear_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Stack_colinear_nc.Update()
c_Stack_colinear_nc.Write()

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
def write_canvas(h1, name, draw_option):
    c_name = "c" + name[1:]
    canvas = rt.TCanvas(str(c_name))
    # h_temp.SetFillColor(rt.kViolet-4)
    h1.Draw(str(draw_option))
    canvas.Update()
    canvas.Draw(str(draw_option))
    canvas.Write()

# write_canvas(h_Reco_colinear_cc, "h_Reco_colinear_cc", "HIST")
# write_canvas(h_misID_colinear_cc, "h_misID_colinear_cc", "HIST")
# write_canvas(h_Sec_colinear_cc, "h_Sec_colinear_cc", "HIST")
write_canvas(h_Prim_colinear_cc, "h_Prim_colinear_cc", "HIST")
# write_canvas(h_Missing_colinear_cc, "h_Missing_colinear_cc", "HIST")
# write_canvas(h_Truecolinear_cc, "h_Truecolinear_cc", "HIST")
write_canvas(h_Effcolinear_cc, "h_Effcolinear_cc", "HIST")

write_canvas(h_Hist_reco_cc, "h_Hist_reco_cc", "HIST")
write_canvas(h_Hist_misID_cc, "h_Hist_misID_cc", "HIST")
write_canvas(h_Hist_sec_cc, "h_Hist_sec_cc", "HIST")
write_canvas(h_Hist_prim_cc, "h_Hist_prim_cc", "HIST")
write_canvas(h_Hist_missing_cc, "h_Hist_missing_cc", "HIST")

write_canvas(h_Uncontained_cc, "h_Uncontained_cc", "HIST")
write_canvas(h_Range_cc, "h_Range_cc", "HIST")
write_canvas(h_Flight_decay_cc, "h_Flight_decay_cc", "HIST")
write_canvas(h_Sec_int_cc, "h_Sec_int_cc", "HIST")


# write_canvas(h_Reco_colinear_nc, "h_Reco_colinear_nc", "HIST")
# write_canvas(h_misID_colinear_nc, "h_misID_colinear_nc", "HIST")
# write_canvas(h_Sec_colinear_nc, "h_Sec_colinear_nc", "HIST")
write_canvas(h_Prim_colinear_nc, "h_prim_colinear_nc", "HIST")
# write_canvas(h_Missing_colinear_nc, "h_Missing_colinear_nc", "HIST")
# write_canvas(h_Truecolinear_nc, "h_Truecolinear_nc", "HIST")
write_canvas(h_Effcolinear_nc, "h_Effcolinear_nc", "HIST")

write_canvas(h_Hist_reco_nc, "h_Hist_reco_nc", "HIST")
write_canvas(h_Hist_misID_nc, "h_Hist_misID_nc", "HIST")
write_canvas(h_Hist_sec_nc, "h_Hist_sec_nc", "HIST")
write_canvas(h_Hist_prim_nc, "h_Hist_prim_nc", "HIST")
write_canvas(h_Hist_missing_nc, "h_Hist_missing_nc", "HIST")

write_canvas(h_Uncontained_nc, "h_Uncontained_nc", "HIST")
write_canvas(h_Range_nc, "h_Range_nc", "HIST")
write_canvas(h_Flight_decay_nc, "h_Flight_decay_nc", "HIST")
write_canvas(h_Sec_int_nc, "h_Sec_int_nc", "HIST")

write_canvas(h_range_detect_eff_cc, "h_range_detect_eff_cc", "HIST")
write_canvas(h_range_detect_eff_nc, "h_range_detect_eff_nc", "HIST")

print(detect_count[0], detect_count[1],detect_count[1]/detect_count[0])
print(prim_counter)
print(missed_long_pions_counter)
f = open("missed_long_track_pi.txt", "w")
f.write(str(missed_long_pions))
f.close()
