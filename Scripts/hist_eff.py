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
        elif (p_Reco_status == 1):
            h_misID_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_cc.Fill(p_History, p_EventWeight)
        elif (p_Reco_status == 2):
            h_Sec_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_cc.Fill(p_History, p_EventWeight)
        elif (p_Reco_status == 3):
            h_Missing_colinear_cc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_cc.Fill(p_History, p_EventWeight)

    if (p_CCNC == 0):
        h_Truecolinear_nc.Fill(p_colinearAng, p_EventWeight)
        if (p_Reco_status == 0):
            h_Reco_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_reco_nc.Fill(p_History, p_EventWeight)
        elif (p_Reco_status == 1):
            h_misID_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_misID_nc.Fill(p_History, p_EventWeight)
        elif (p_Reco_status == 2):
            h_Sec_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_sec_nc.Fill(p_History, p_EventWeight)
        elif (p_Reco_status == 3):
            h_Missing_colinear_nc.Fill(p_colinearAng, p_EventWeight)
            h_Hist_missing_nc.Fill(p_History, p_EventWeight)


h_Reco_colinear_cc.Scale(6.67e20/t_pot)
h_misID_colinear_cc.Scale(6.67e20/t_pot)
h_Sec_colinear_cc.Scale(6.67e20/t_pot)
h_Missing_colinear_cc.Scale(6.67e20/t_pot)
h_Truecolinear_cc.Scale(6.67e20/t_pot)
h_Hist_reco_cc.Scale(6.67e20/t_pot)
h_Hist_misID_nc.Scale(6.67e20/t_pot)
h_Hist_sec_cc.Scale(6.67e20/t_pot)
h_Hist_missing_cc.Scale(6.67e20/t_pot)
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
h_Effcolinear_nc.Divide(h_Reco_colinear_nc, h_Truecolinear_nc)

c_Stack_colinear_cc = rt.TCanvas("c_Stack_colinear_cc")
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Reco_colinear_cc
h_temp.SetFillColor(rt.kGreen-3)
h_temp.Draw("bar")
h_Stack_colinear_cc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_misID_colinear_cc
h_temp.SetFillColor(rt.kBlue-7)
h_temp.Draw("bar")
h_Stack_colinear_cc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Sec_colinear_cc
h_temp.SetFillColor(rt.kYellow-7)
h_temp.Draw("bar")
h_Stack_colinear_cc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Missing_colinear_cc
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.1, 0.55, 0.3, 0.9)
legend.AddEntry(h_Reco_colinear_cc,"Good Reco","f")
legend.AddEntry(h_misID_colinear_cc,"Wrong PDG","f")
legend.AddEntry(h_Sec_colinear_cc,"Secondary Track","f")
legend.AddEntry(h_temp,"Missing Track","f")
h_temp_leg = h_temp.Clone("h_temp_leg")
h_temp_leg.GetListOfFunctions().Add(legend)
h_temp_leg.Draw("bar")
legend.Draw()
h_Stack_colinear_cc.Add(h_temp_leg)
h_Stack_colinear_cc.Draw("bar")
c_Stack_colinear_cc.Update()
c_Stack_colinear_cc.Write()

c_Stack_colinear_nc = rt.TCanvas("c_Stack_colinear_nc")
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Reco_colinear_nc
h_temp.SetFillColor(rt.kGreen-3)
h_temp.Draw("bar")
h_Stack_colinear_nc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_misID_colinear_nc
h_temp.SetFillColor(rt.kBlue-7)
h_temp.Draw("bar")
h_Stack_colinear_nc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Sec_colinear_nc
h_temp.SetFillColor(rt.kYellow-7)
h_temp.Draw("bar")
h_Stack_colinear_nc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1)
h_temp = h_Missing_colinear_nc
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.1, 0.55, 0.3, 0.9)
legend.AddEntry(h_Reco_colinear_nc,"Good Reco","f")
legend.AddEntry(h_misID_colinear_nc,"Wrong PDG","f")
legend.AddEntry(h_Sec_colinear_nc,"Secondary Track","f")
legend.AddEntry(h_temp,"Missing Track","f")
h_temp_leg = h_temp.Clone("h_temp_leg")
h_temp_leg.GetListOfFunctions().Add(legend)
h_temp_leg.Draw("bar")
legend.Draw()
h_Stack_colinear_nc.Add(h_temp_leg)
h_Stack_colinear_nc.Draw("bar")
c_Stack_colinear_nc.Update()
c_Stack_colinear_nc.Write()

h_Reco_colinear_cc.Draw()
h_misID_colinear_cc.Draw()
h_Sec_colinear_cc.Draw()
h_Missing_colinear_cc.Draw()
h_Truecolinear_cc.Draw()
h_Effcolinear_cc.Draw()
h_Hist_reco_cc.Draw()
h_Hist_misID_nc.Draw()
h_Hist_sec_cc.Draw()
h_Hist_missing_cc.Draw()

h_Reco_colinear_nc.Draw()
h_misID_colinear_nc.Draw()
h_Sec_colinear_nc.Draw()
h_Missing_colinear_nc.Draw()
h_Truecolinear_nc.Draw()
h_Effcolinear_nc.Draw()
h_Hist_reco_nc.Draw()
h_Hist_misID_nc.Draw()
h_Hist_sec_nc.Draw()
h_Hist_missing_nc.Draw()

f_out.cd()
h_Reco_colinear_cc.Write()
h_misID_colinear_cc.Write()
h_Sec_colinear_cc.Write()
h_Missing_colinear_cc.Write()
h_Truecolinear_cc.Write()
h_Effcolinear_cc.Write()
h_Hist_reco_cc.Write()
h_Hist_misID_nc.Write()
h_Hist_sec_cc.Write()
h_Hist_missing_cc.Write()
h_Stack_colinear_cc.Write()

h_Reco_colinear_nc.Write()
h_misID_colinear_nc.Write()
h_Sec_colinear_nc.Write()
h_Missing_colinear_nc.Write()
h_Truecolinear_nc.Write()
h_Effcolinear_nc.Write()
h_Hist_reco_nc.Write()
h_Hist_misID_nc.Write()
h_Hist_sec_nc.Write()
h_Hist_missing_nc.Write()
h_Stack_colinear_nc.Write()
f_out.Close()


