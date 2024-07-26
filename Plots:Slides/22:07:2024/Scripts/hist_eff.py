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

t_pot = 0.
for entry in p1:
    t_pot += entry.totGoodPOT

h_Reco_colinear_cc = rt.TH1F("h_Reco_colinear_cc","Number of Reconstructed Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_misID_colinear_cc = rt.TH1F("h_misID_colinear_cc","Number of mis-ID'd Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Sec_colinear_cc = rt.TH1F("h_Sec_colinear_cc","Number Pions Reco'd as Secondary by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Missing_colinear_cc = rt.TH1F("h_Missing_colinear_cc","Number of Missing Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Truecolinear_cc = rt.TH1F("h_Truecolinear_cc","Number of True Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Effcolinear_cc = rt.TH1F("h_Effcolinear_cc","Reco Efficiency of Pions by Most Colinear Primary (cc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Stack_colinear_cc = rt.THStack("h_Stack_colinear_cc", "Reco Status of Pions by Most Colinear Primary (cc)")

h_Hist_reco_cc = rt.TH1F("h_Hist_reco_cc", "History of Reconstructed Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_misID_cc = rt.TH1F("h_Hist_misID_cc", "History of misID'd Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_sec_cc = rt.TH1F("h_Hist_sec_cc", "History of Pions ID'd as Secondary (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_missing_cc = rt.TH1F("h_Hist_missing_cc", "History of Missing Pions (cc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
# h_Hist_stack_cc = rt.THStack("h_Hist_stack_cc", "History of Pions by Reco Type (cc)")

h_Uncontained_cc = rt.TH1F("h_Uncontained_cc", "Reco Result of Uncontained Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Range_cc = rt.TH1F("h_Range_cc", "Reco Result of Ranged out Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Flight_decay_cc = rt.TH1F("h_Flight_decay_cc", "Reco Result of Decayed in Flight Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Sec_int_cc = rt.TH1F("h_Sec_int_cc", "Reco Result of Secondary Interacted Pions (cc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)

h_Dist_reco_cc = rt.TH1F("h_Dist_reco_cc", "Distance to Interaction of Reconstructed Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_misID_cc = rt.TH1F("h_Dist_misID_cc", "Distance to Interaction of mis-ID Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_sec_cc = rt.TH1F("h_Dist_sec_cc", "Distance to Interaction of ID'd as Sec. Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_missing_cc = rt.TH1F("h_Dist_missing_cc", "Distance to Interaction of Missing Pions (cc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_stack_cc = rt.THStack("h_Dist_stack_cc", "Distance to Interaction of Pions by Reco Type (cc)")

h_Dist_reco_hist_cc = rt.THStack("h_Dist_reco_hist_cc", "Distacce to Interaction of Reconstructed Pions by History (cc)")
h_Dist_misID_hist_cc = rt.THStack("h_Dist_misID_hist_cc", "Distacce to Interaction of mis-ID Pions by History (cc)")
h_Dist_sec_hist_cc = rt.THStack("h_Dist_sec_hist_cc", "Distacce to Interaction of ID'd as Sec. Pions by History (cc)")
h_Dist_missing_hist_cc = rt.THStack("h_Dist_missing_hist_cc", "Distacce to Interaction of Missing Pions by History (cc)")
Dist_hist_reco_list_cc = [h_Dist_reco_hist_cc,h_Dist_misID_hist_cc,h_Dist_sec_hist_cc,h_Dist_missing_hist_cc]

h_Reco_colinear_nc = rt.TH1F("h_Reco_colinear_nc","Number of Reconstructed Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_misID_colinear_nc = rt.TH1F("h_misID_colinear_nc","Number of mis-ID'd Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Sec_colinear_nc = rt.TH1F("h_Sec_colinear_nc","Number Pions Reco'd as Secondary by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Missing_colinear_nc = rt.TH1F("h_Missing_colinear_nc","Number of Missing Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Truecolinear_nc = rt.TH1F("h_Truecolinear_nc","Number of True Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Effcolinear_nc = rt.TH1F("h_Effcolinear_nc","Reco Efficiency of Pions by Most Colinear Primary (nc);Colinear Primary Angle (cos theta);# of Pions per 6.67e20 POT", 40,-1,1)
h_Stack_colinear_nc = rt.THStack("h_Stack_colinear_nc", "Reco Status of Pions by Most Colinear Primary (nc)")

h_Hist_reco_nc = rt.TH1F("h_Hist_reco_nc", "History of Reconstructed Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_misID_nc = rt.TH1F("h_Hist_misID_nc", "History of misID'd Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_sec_nc = rt.TH1F("h_Hist_sec_nc", "History of Pions ID'd as Secondary (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
h_Hist_missing_nc = rt.TH1F("h_Hist_missing_nc", "History of Missing Pions (nc);Pion History Number;# of Pions per 6.67e20 POT", 4,0,4)
# h_Hist_stack_nc = rt.THStack("h_Hist_stack_nc", "History of Pions by Reco Type (cc)")

h_Uncontained_nc = rt.TH1F("h_Uncontained_nc", "Reco Result of Uncontained Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Range_nc = rt.TH1F("h_Range_nc", "Reco Result of Ranged out Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Flight_decay_nc = rt.TH1F("h_Flight_decay_nc", "Reco Result of Decayed in Flight Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)
h_Sec_int_nc = rt.TH1F("h_Sec_int_nc", "Reco Result of Secondary Interacted Pions (nc);Reco Type; # of Pions per 6.67e20 POT", 4,0,4)

h_Dist_reco_nc = rt.TH1F("h_Dist_reco_nc", "Distance to Interaction of Reconstructed Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_misID_nc = rt.TH1F("h_Dist_misID_nc", "Distance to Interaction of mis-ID Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_sec_nc = rt.TH1F("h_Dist_sec_nc", "Distance to Interaction of ID'd as Sec. Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_missing_nc = rt.TH1F("h_Dist_missing_nc", "Distance to Interaction of Missing Pions (nc);Distance to Interaction (cm);# of Pions per 6.67e20 POT", 100,0,50)
h_Dist_stack_nc = rt.THStack("h_Dist_stack_nc", "Distance to Interaction of Pions by Reco Type (nc)")

h_Dist_reco_hist_nc = rt.THStack("h_Dist_reco_hist_nc", "Distance to Interaction of Reconstructed Pions by History (nc)")
h_Dist_misID_hist_nc = rt.THStack("h_Dist_misID_hist_nc", "Distance to Interaction of mis-ID Pions by History (nc)")
h_Dist_sec_hist_nc = rt.THStack("h_Dist_sec_hist_nc", "Distance to Interaction of ID'd as Sec. Pions by History (nc)")
h_Dist_missing_hist_nc = rt.THStack("h_Dist_missing_hist_nc", "Distance to Interaction of Missing Pions by History (nc)")
Dist_hist_reco_list_nc = [h_Dist_reco_hist_nc,h_Dist_misID_hist_nc,h_Dist_sec_hist_nc,h_Dist_missing_hist_nc]
temp_hists = []
for i in range(0,32):
    h_temp = rt.TH1F("","",100,0,100)
    temp_hists.append(h_temp)

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

    # CC pions
    if (p_CCNC == 0):
        h_Truecolinear_cc.Fill(p_colinearAng, p_EventWeight)
        if (p_Reco_status == 0):
            h_Reco_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_cc.Fill(p_History, p_EventWeight)
            h_Dist_reco_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 1):
            h_misID_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_cc.Fill(p_History, p_EventWeight)
            h_Dist_misID_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 2):
            h_Sec_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_cc.Fill(p_History, p_EventWeight)
            h_Dist_sec_cc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 3):
            h_Missing_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_cc.Fill(p_History, p_EventWeight)
            h_Dist_missing_cc.Fill(p_Dist_to_interaction, p_EventWeight)
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
            h_Uncontained_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 1):
            h_Range_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 2):
            h_Flight_decay_cc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 3):
            h_Sec_int_cc.Fill(p_Reco_status, p_EventWeight)
    # NC pions
    if (p_CCNC == 1):
        h_Truecolinear_nc.Fill(p_colinearAng, p_EventWeight)
        if (p_Reco_status == 0):
            h_Reco_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_nc.Fill(p_History, p_EventWeight)
            h_Dist_reco_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status + 16
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 1):
            h_misID_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_nc.Fill(p_History, p_EventWeight)
            h_Dist_misID_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status + 16
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 2):
            h_Sec_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_nc.Fill(p_History, p_EventWeight)
            h_Dist_sec_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status + 16
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        elif (p_Reco_status == 3):
            h_Missing_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_nc.Fill(p_History, p_EventWeight)
            h_Dist_missing_nc.Fill(p_Dist_to_interaction, p_EventWeight)
            cat = 4*p_Reco_status + 16
            if (p_History == 0):
                temp_hists[cat].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 1):
                temp_hists[cat + 1].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 2):
                temp_hists[cat + 2].Fill(p_Dist_to_interaction, p_EventWeight)
            elif (p_History == 3):
                temp_hists[cat + 3].Fill(p_Dist_to_interaction, p_EventWeight)
        if (p_History == 0):
            h_Uncontained_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 1):
            h_Range_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 2):
            h_Flight_decay_nc.Fill(p_Reco_status, p_EventWeight)
        elif (p_History == 3):
            h_Sec_int_nc.Fill(p_Reco_status, p_EventWeight)

h_Reco_colinear_cc.Scale(6.67e20/t_pot)
h_misID_colinear_cc.Scale(6.67e20/t_pot)
h_Sec_colinear_cc.Scale(6.67e20/t_pot)
h_Missing_colinear_cc.Scale(6.67e20/t_pot)
h_Truecolinear_cc.Scale(6.67e20/t_pot)
h_Hist_reco_cc.Scale(6.67e20/t_pot)
h_Hist_misID_cc.Scale(6.67e20/t_pot)
h_Hist_sec_cc.Scale(6.67e20/t_pot)
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

for i in range(0,32):
    temp_hists[i].Scale(6.67e20/t_pot)
    print(temp_hists[i].GetEntries())

History_tags = ["Uncontained", "Range Out", "Sec. Interaction", "Decay in Flight"]
Reco_tags = ["Good Reco", "Mis-ID", "ID'd as Sec.", "Missing"]
def write_stack(h1, h2, h3, h4, h_stack, name, l_pos, l_tags):
    canvas = rt.TCanvas(str(name))
    h_temp = h1
    h_temp.SetFillColor(rt.kGreen-3)
    h_temp.Draw("HIST")
    h_stack.Add(h_temp)
    h_temp = h2
    h_temp.SetFillColor(rt.kBlue-7)
    h_temp.Draw("HIST")
    h_stack.Add(h_temp)
    h_temp = h3
    h_temp.SetFillColor(rt.kYellow-7)
    h_temp.Draw("HIST")
    h_stack.Add(h_temp)
    h_temp = h4
    h_temp.SetFillColor(rt.kRed-7)
    legend = rt.TLegend(l_pos[0], l_pos[1], l_pos[2], l_pos[3])
    legend.AddEntry(h1,l_tags[0],"f")
    legend.AddEntry(h2,l_tags[1],"f")
    legend.AddEntry(h3,l_tags[2],"f")
    legend.AddEntry(h_temp,l_tags[3],"f")
    h_temp_leg = h_temp.Clone("h_temp_leg")
    h_temp_leg.GetListOfFunctions().Add(legend)
    h_temp_leg.Draw()
    legend.Draw()
    h_stack.Add(h_temp_leg)
    h_stack.Draw("HIST")
    return canvas
f_out.cd()
c_Stack_colinear_cc = write_stack(h_Reco_colinear_cc,h_misID_colinear_cc,
                                  h_Sec_colinear_cc,h_Missing_colinear_cc,
                                  h_Stack_colinear_cc,"c_Stack_colinear_cc",
                                  [0.1, 0.55, 0.3, 0.9], Reco_tags)
h_Stack_colinear_cc.GetXaxis().SetTitle("Most Colinear Primary Angle (cos theta)")
h_Stack_colinear_cc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Stack_colinear_cc.Update()
c_Stack_colinear_cc.Write()

c_Dist_stack_cc = write_stack(h_Dist_reco_cc,h_Dist_misID_cc,h_Dist_sec_cc,
                              h_Dist_missing_cc,h_Dist_stack_cc,"c_Dist_stack_cc",
                              [0.7, 0.55, 0.9, 0.9], Reco_tags)
h_Dist_stack_cc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_cc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_cc.Update()
c_Dist_stack_cc.Write()

c_Stack_colinear_nc = write_stack(h_Reco_colinear_nc,h_misID_colinear_nc,
                                  h_Sec_colinear_nc,h_Missing_colinear_nc,
                                  h_Stack_colinear_nc,"c_Stack_colinear_nc",
                                  [0.1, 0.55, 0.3, 0.9], Reco_tags)
h_Stack_colinear_nc.GetXaxis().SetTitle("Most Colinear Primary Angle (cos theta)")
h_Stack_colinear_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Stack_colinear_nc.Update()
c_Stack_colinear_nc.Write()

c_Dist_stack_nc = write_stack(h_Dist_reco_nc,h_Dist_misID_nc,h_Dist_sec_nc,
                              h_Dist_missing_nc,h_Dist_stack_nc,"c_Dist_stack_nc",
                              [0.7, 0.55, 0.9, 0.9], Reco_tags)
h_Dist_stack_nc.GetXaxis().SetTitle("Distance to Interaction (cm)")
h_Dist_stack_nc.GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
c_Dist_stack_nc.Update()
c_Dist_stack_nc.Write()

nc_stack_names = ["c_Dist_reco_hist_nc","c_Dist_misID_hist_nc","c_Dist_sec_hist_nc","c_Dist_missing_hist_nc"]
for i in range(4,8):
    canvas = write_stack(temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3],
                Dist_hist_reco_list_nc[i-4],nc_stack_names[i-4],[0.7, 0.55, 0.9, 0.9], History_tags)
    Dist_hist_reco_list_nc[i-4].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_nc[i-4].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

cc_stack_names = ["c_Dist_reco_hist_cc","c_Dist_misID_hist_cc","c_Dist_sec_hist_cc","c_Dist_missing_hist_cc"]
for i in range(0,4):
    canvas = write_stack(temp_hists[4*i], temp_hists[4*i+1], temp_hists[4*i+2],temp_hists[4*i+3],
                Dist_hist_reco_list_cc[i],cc_stack_names[i],[0.7, 0.55, 0.9, 0.9], History_tags)
    Dist_hist_reco_list_cc[i].GetXaxis().SetTitle("Distance to Interaction (cm)")
    Dist_hist_reco_list_cc[i].GetYaxis().SetTitle("# of Pions per 6.67e20 POT")
    canvas.Update()
    canvas.Write()

h_Reco_colinear_cc.Draw()
h_misID_colinear_cc.Draw()
h_Sec_colinear_cc.Draw()
h_Missing_colinear_cc.Draw()
h_Truecolinear_cc.Draw()
h_Effcolinear_cc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_reco_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_reco_cc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_misID_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_misID_cc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_sec_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_sec_cc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_missing_cc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_missing_cc.Draw()

for i in range(0, len(History_tags)):
    h_Uncontained_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Uncontained_cc.Draw()
for i in range(0, len(History_tags)):
    h_Range_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Range_cc.Draw()
for i in range(0, len(History_tags)):
    h_Flight_decay_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Flight_decay_cc.Draw()
for i in range(0, len(History_tags)):
    h_Sec_int_cc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Sec_int_cc.Draw()

h_Reco_colinear_nc.Draw()
h_misID_colinear_nc.Draw()
h_Sec_colinear_nc.Draw()
h_Missing_colinear_nc.Draw()
h_Truecolinear_nc.Draw()
h_Effcolinear_nc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_reco_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_reco_nc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_misID_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_misID_nc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_sec_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_sec_nc.Draw()
for i in range(0, len(History_tags)):
    h_Hist_missing_nc.GetXaxis().SetBinLabel(i+1, History_tags[i])
h_Hist_missing_nc.Draw()

for i in range(0, len(History_tags)):
    h_Uncontained_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Uncontained_nc.Draw()
for i in range(0, len(History_tags)):
    h_Range_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Range_nc.Draw()
for i in range(0, len(History_tags)):
    h_Flight_decay_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Flight_decay_nc.Draw()
for i in range(0, len(History_tags)):
    h_Sec_int_nc.GetXaxis().SetBinLabel(i+1, Reco_tags[i])
h_Sec_int_nc.Draw()

def write_canvas(h1, name):
    c_name = "c" + name[1:]
    canvas = rt.TCanvas(str(c_name))
    h1.Draw("bar")
    canvas.Update()
    canvas.Draw("bar")
    canvas.Write()
    
# write_canvas(h_Reco_colinear_cc, "h_Reco_colinear_cc")
# write_canvas(h_misID_colinear_cc, "h_misID_colinear_cc")
# write_canvas(h_Sec_colinear_cc, "h_Sec_colinear_cc")
# write_canvas(h_Missing_colinear_cc, "h_Missing_colinear_cc")
# write_canvas(h_Truecolinear_cc, "h_Truecolinear_cc")
write_canvas(h_Effcolinear_cc, "h_Effcolinear_cc")

write_canvas(h_Hist_reco_cc, "h_Hist_reco_cc")
write_canvas(h_Hist_misID_cc, "h_Hist_misID_cc")
write_canvas(h_Hist_sec_cc, "h_Hist_sec_cc")
write_canvas(h_Hist_missing_cc, "h_Hist_missing_cc")

write_canvas(h_Uncontained_cc, "h_Uncontained_cc")
write_canvas(h_Range_cc, "h_Range_cc")
write_canvas(h_Flight_decay_cc, "h_Flight_decay_cc")
write_canvas(h_Sec_int_cc, "h_Sec_int_cc")


# write_canvas(h_Reco_colinear_nc, "h_Reco_colinear_nc")
# write_canvas(h_misID_colinear_nc, "h_misID_colinear_nc")
# write_canvas(h_Sec_colinear_nc, "h_Sec_colinear_nc")
# write_canvas(h_Missing_colinear_nc, "h_Missing_colinear_nc")
# write_canvas(h_Truecolinear_nc, "h_Truecolinear_nc")
write_canvas(h_Effcolinear_nc, "h_Effcolinear_nc")

write_canvas(h_Hist_reco_nc, "h_Hist_reco_nc")
write_canvas(h_Hist_misID_nc, "h_Hist_misID_nc")
write_canvas(h_Hist_sec_nc, "h_Hist_sec_nc")
write_canvas(h_Hist_missing_nc, "h_Hist_missing_nc")

write_canvas(h_Uncontained_nc, "h_Uncontained_nc")
write_canvas(h_Range_nc, "h_Range_nc")
write_canvas(h_Flight_decay_nc, "h_Flight_decay_nc")
write_canvas(h_Sec_int_nc, "h_Sec_int_nc")
f_out.Close()


