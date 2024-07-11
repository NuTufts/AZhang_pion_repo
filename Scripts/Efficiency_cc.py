# Efficiency_cc.py
# Author: Andy Zhang
# Purpose: Find Reconstruction Efficiencies of True Primary Pions in CC Events ordered by Kinetic Energy, categorized by other Primary particles. Also organize improperly reconstructed pions by cause (Incorrect PDG, missing track, reconstructed as secondary).
# Data Created: June 5th, 2024
# Last Modified: 3rd, 2024

import ROOT as rt
import numpy as np
from ROOT import gROOT, AddressOf
f1 = rt.TFile("inputFile.root")
t1 = f1.Get("EventTree")
p1 = f1.Get("potTree")

class pion_container:
    def __init__(self, true_trackID, top_recoComp, true_pionE, xsecWeight):
        self.true_trackID = true_trackID
        self.reco_found = False
        self.TID_found = False
        self.top_recoComp = top_recoComp
        self.true_pionE = true_pionE
        self.xsecWeight = xsecWeight
    def toggle_reco_found(self):
        self.reco_found = True
    def toggle_TID_found(self):
        self.TID_found = True
    def set_top_recoComp(self, trackTrueComp):
        self.top_recoComp = trackTrueComp
        
# Recursively Test
def prim_ancestor(trueSimPartTID, target_TID):
    index = 999
    for i in range(0,len(trueSimPartTID)):
        if (trueSimPartTID[i] == target_TID):
            index = i
    if(index == 999):
#        print("Fatal.")
        return index
    if(trueSimPartProcess[index] == 0):
#        print("Primary.")
        return index
        
    else:
#        print("Secondary...")
        return prim_ancestor(trueSimPartTID, trueSimPartMID[index])
        

# Count total POT
t_pot = 0.
for entry in p1:
    t_pot += entry.totGoodPOT

h_simE_cc = rt.TH1F("h_simE_cc","Number of Reconstructed Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_trueE_cc = rt.TH1F("h_trueE_cc","Number of Truth-Matched Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_Eff_cc = rt.TH1F("h_Eff_cc","Pion Reconstruction Efficiency by Energy (cc);Pion Kinetic Energy (MeV);Reconstruction Efficiency", 40,0,1200)

h_recoComp_cc = rt.TH1F("h_recoComp_cc","Greatest Track Completeness of Truth-Matched Pion Events (cc);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

h_simE_Mu = rt.TH1F("h_simE_Mu","Number of Reconstructed Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_trueE_Mu = rt.TH1F("h_trueE_Mu","Number of Truth-Matched Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_Eff_Mu = rt.TH1F("h_Eff_Mu","Pion Reconstruction Efficiency by Energy with Muons (cc);Pion Kinetic Energy (MeV);Reconstruction Efficiency", 40,0,1200)

h_recoComp_Mu = rt.TH1F("h_recoComp_Mu","Greatest Track Completeness of Truth-Matched Pion Events (cc);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

h_simE_Pr = rt.TH1F("h_simE_Pr","Number of Reconstructed Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_trueE_Pr = rt.TH1F("h_trueE_Pr","Number of Truth-Matched Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_Eff_Pr = rt.TH1F("h_Eff_Pr","Pion Reconstruction Efficiency by Energy with Protons (cc);Pion Kinetic Energy (MeV);Reconstruction Efficiency", 40,0,1200)

h_recoComp_Pr = rt.TH1F("h_recoComp_Pr","Greatest Track Completeness of Truth-Matched Pion Events (cc);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

h_simE_Neu = rt.TH1F("h_simE_Neu","Number of Reconstructed Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_trueE_Neu = rt.TH1F("h_trueE_Neu","Number of Truth-Matched Pions by Energy (cc);Pion Kinetic Energy (MeV);# of Pions per 6.67e20 POT", 40,0,1200)

h_Eff_Neu = rt.TH1F("h_Eff_Neu","Pion Reconstruction Efficiency by Energy with Neutrons (cc);Pion Kinetic Energy (MeV);Reconstruction Efficiency", 40,0,1200)

h_recoComp_Neu = rt.TH1F("h_recoComp_Neu","Greatest Track Completeness of Truth-Matched Pion Events (cc);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

h_missing = rt.TH1F("h_missing","Improperly Reconstructed Pions with Missing Tracks by Energy (cc);Pion Kinetic Energy (MeV);# of Mis-Identified Tracks", 40,0,1200)

h_secondary = rt.TH1F("h_secondary","Improperly Reconstructed Pions Mis-Identified as their Secondaries by Energy (cc);Pion Kinetic Energy (MeV);# of Mis-Identified Tracks", 40,0,1200)

h_misID = rt.TH1F("h_misID","Improperly Reconstructed Pion Mis-Identified as Other Primary Particles by Energy (cc);Pion Kinetic Energy (MeV);# of Improperly Reconstructed Pions", 40,0,1200)

h_improper = rt.TH1F("h_improper","Total Improperly Reconstructed Pions by Energy (cc); Pion Kinetic Energy (MeV); # of Improperly Reconstructed Pions", 40,0,1200)

h_unconstructed_stack = rt.THStack("h_unconstructed_stack","Cause of Improper Pion Reconstruction by Energy (cc)")

cc = [] # List for event storage
# A bunch of checks and counters
entry_count = 0 # Number of entries
pion_count = [0,0,0] # Number of True Pions; True Pions with a reco track of same TID; True Pions with a reco track of same TID and PDG
errors = 0 # Number of Entries where the number of pions without a reco track != missing or mis-ID'd-as-secondary pions (due to duplicate tracks)
filled_secondary = 0 # Fill counter for the Mis-ID'd-as-Secondary histogram
filled_missing = 0 # Fill counter for the missing pion histogram
dupe_counter = 0 # Number of duplicate tracks for the same pion across all entries
total, toggle = 0,0 # Total number of pions in file; total number of reconstructed pion tracks in file
filtered = [0,0,0,0,0,0] # Number of pions whose interaction involved primary Muons, Protons, Neutrons; number of said pions properly reconstructed

for entry in t1:
    xsecWeight = entry.xsecWeight
    trueNuPDG = entry.trueNuPDG
    nTracks = entry.nTracks
    trackTrueTID = entry.trackTrueTID
    nTrueSimParts = entry.nTrueSimParts
    trueSimPartTID = entry.trueSimPartTID
    trueSimPartMID = entry.trueSimPartMID
    trueSimPartE = entry.trueSimPartE
    trueSimPartPx = entry.trueSimPartPx
    trueSimPartPy = entry.trueSimPartPy
    trueSimPartPz = entry.trueSimPartPz
    trueSimPartPDG = entry.trueSimPartPDG
    trueSimPartProcess = entry.trueSimPartProcess
    trueNuCCNC = entry.trueNuCCNC
    trackPID = entry.trackPID
    trackTrueComp = entry.trackTrueComp
    trackIsSecondary = entry.trackIsSecondary

    event_secondary, event_missing = 0, 0
    # Filter for CC Events
    if(trueNuCCNC == 0):
        event_pions = [] # List for pions in event
        event_other = [0, 0, 0] # List for other primaries in event
        event_pion_count = [0,0] # Same as pion_count, but for each event
        
        # Create container class for each pion in event
        for j in range(0, nTrueSimParts):
            KE = trueSimPartE[j] - np.sqrt(trueSimPartE[j]**2 - (trueSimPartPx[j]**2 + trueSimPartPy[j]**2 + trueSimPartPz[j]**2))
            if(trueSimPartProcess[j] == 0):
                if(trueSimPartPDG[j] == 211 or trueSimPartPDG[j] == -211):
                    if(KE >= 30):
                        temp = pion_container(trueSimPartTID[j], 0, KE, xsecWeight)
                        event_pions.append(temp)
                # Count other primary particles
                if(trueSimPartPDG[j] == 13 and KE >= 30):
                    event_other[0] += 1
                if(trueSimPartPDG[j] == 2212 and KE >= 60):
                    event_other[1] += 1
                if(trueSimPartPDG[j] == 2112 and KE >= 60):
                    event_other[2] += 1
        
        # Find reconstructed pions by matching track IDs. Also find maximum track completeness
        for j in range(0, len(event_pions)):
            pion_count[0] += 1
            event_pion_count[0] += 1
            primary_index = []
            for k in range(0, nTracks):
                if(trackIsSecondary[k] != 1):
                    if(event_pions[j].true_trackID == trackTrueTID[k]):
                        # If triggered, a track with the same TID was already found, and there is a duplicate
                        if(event_pions[j].TID_found == True):
                            dupe_counter += 1
                        event_pions[j].toggle_TID_found()
                        event_pion_count[1] += 1
                        if(trackPID[k] == 211 or trackPID[k] == -211):
                            event_pions[j].toggle_reco_found()
                            # Fill the highest Completeness
                            if(trackTrueComp[k] > event_pions[j].top_recoComp):
                                event_pions[j].set_top_recoComp(trackTrueComp[k])
                
                # Recursively find ancestor particle's index
                track_ancestor = prim_ancestor(trueSimPartTID, trackTrueTID[k])
                primary_index.append(track_ancestor)
                
            # Fill histograms according to pion recostruction type
            if(event_pions[j].TID_found == True):
                pion_count[1] += 1
                if(event_pions[j].reco_found == False): # Track was found, but incorrect PDG
                    h_misID.Fill(event_pions[j].true_pionE, event_pions[j].xsecWeight)
                else:
                    pion_count[2] += 1
            else:
                primary_found = False
                for i in range(0, len(primary_index)):
                    if(primary_index[i] == 999):
                        continue
                    # Track not found, but a secondary of the pion was found
                    if(trueSimPartTID[primary_index[i]] == event_pions[j].true_trackID):
                        h_secondary.Fill(event_pions[j].true_pionE, event_pions[j].xsecWeight)
                        primary_found = True
                        filled_secondary += 1
                        event_secondary += 1
                        break
                # Track not found, secondaries not found
                if(primary_found == False):
                    h_missing.Fill(event_pions[j].true_pionE, event_pions[j].xsecWeight)
                    event_missing += 1
                    filled_missing += 1
#        if(event_pion_count[0] - event_pion_count[1] != event_secondary + event_missing):
#            errors += 1
#        total += event_pion_count[0]
#        toggle += event_pion_count[1]
        event_pions.append(event_other)
        cc.append(event_pions)
    entry_count += 1
print("Number of Pions with Missing Tracks: ", filled_missing, " Number of Pions identified as their Secondaries: ", filled_secondary)
print("Number of Events where sum of missing and secondaried pions != total trackless pions:", errors, " Number of Duplicate Tracks", dupe_counter)
#print(total, toggle)

# Fill histograms by category
reco_fill_count = 0 # Counter for filled reconstructed pions
for i in range(0, len(cc)):
    pion_index = len(cc[i]) - 1
    for j in range(0, pion_index):
        h_trueE_cc.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
        if(cc[i][pion_index][0] >= 1):
            h_trueE_Mu.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
            filtered[0] += 1
        if(cc[i][pion_index][1] >= 1):
            h_trueE_Pr.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
            filtered[1] += 1
        if(cc[i][pion_index][2] >= 1):
            h_trueE_Neu.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
            filtered[2] += 1
        if(cc[i][j].reco_found == True):
            h_recoComp_cc.Fill(cc[i][j].top_recoComp)
            h_simE_cc.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
            reco_fill_count += 1
            if(cc[i][pion_index][0] >= 1):
                h_simE_Mu.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
                filtered[3] += 1
                h_recoComp_Mu.Fill(cc[i][j].top_recoComp)
            if(cc[i][pion_index][1] >= 1):
                h_simE_Pr.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
                filtered[4] += 1
                h_recoComp_Pr.Fill(cc[i][j].top_recoComp)
            if(cc[i][pion_index][2] >= 1):
                h_simE_Neu.Fill(cc[i][j].true_pionE, cc[i][j].xsecWeight)
                filtered[5] += 1
                h_recoComp_Neu.Fill(cc[i][j].top_recoComp)

print("number of entries: ", entry_count)
print("total primary pions: ", pion_count[0])
print("Reco'd Pions: ", pion_count[2])
print("Non-Reco Pions: ")
print("Total -> ", pion_count[0] - pion_count[2], "Track Mis-ID -> ", pion_count[1] - pion_count[2], " or ", pion_count[1] - reco_fill_count, "Missing Track -> ", pion_count[0] - pion_count[1])
print(filtered)

# Scale all histograms to target POT
h_simE_cc.Scale(6.67e20/t_pot)
h_trueE_cc.Scale(6.67e20/t_pot)
h_simE_Mu.Scale(6.67e20/t_pot)
h_trueE_Mu.Scale(6.67e20/t_pot)
h_simE_Pr.Scale(6.67e20/t_pot)
h_trueE_Pr.Scale(6.67e20/t_pot)
h_simE_Neu.Scale(6.67e20/t_pot)
h_trueE_Neu.Scale(6.67e20/t_pot)
h_missing.Scale(6.67e20/t_pot)
h_misID.Scale(6.67e20/t_pot)
h_secondary.Scale(6.67e20/t_pot)

# Create stacks for efficiency
h_stack_cc = rt.THStack("h_stack_cc","Reconstruction of True Pions by Energy (cc)")
h_stack_Mu = rt.THStack("h_stack_Mu","Reconstruction of True Pions with Muons by Energy (cc)")
h_stack_Pr = rt.THStack("h_stack_Pr","Reconstruction of True Pions with Protons by Energy (cc)")
h_stack_Neu = rt.THStack("h_stack_Neu","Reconstruction of True Pions with Neutrons by Energy (cc)")

# Plot Beautification and Writing
f_out = rt.TFile("Pi_Efficiency_cc.root","RECREATE")
c_stack_cc = rt.TCanvas("c_stack_cc")
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_simE_cc
h_temp.SetFillColor(rt.kGreen-3)
h_stack_cc.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp.Add(h_trueE_cc, h_simE_cc, 1, -1)
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.7, 0.1, 0.9, 0.3)
legend.AddEntry(h_simE_cc,"Good Reco","f")
legend.AddEntry(h_temp,"Improper Reco","f")
h_temp.GetListOfFunctions().Add(legend)
h_temp.Draw("bar")
h_stack_cc.Add(h_temp)
legend.Draw()
h_stack_cc.Draw("bar")
c_stack_cc.Update()
c_stack_cc.Write()

c_stack_Mu = rt.TCanvas("c_stack_Mu")
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_simE_Mu
h_temp.SetFillColor(rt.kGreen-3)
h_temp.Draw("bar")
h_stack_Mu.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp.Add(h_trueE_Mu, h_simE_Mu, 1, -1)
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.7, 0.1, 0.9, 0.3)
legend.AddEntry(h_simE_Mu,"Good Reco","f")
legend.AddEntry(h_temp,"Improper Reco","f")
h_temp.GetListOfFunctions().Add(legend)
h_temp.Draw("bar")
h_stack_Mu.Add(h_temp)
legend.Draw()
h_stack_Mu.Draw("bar")
c_stack_Mu.Update()
c_stack_Mu.Write()

c_stack_Pr = rt.TCanvas("c_stack_Pr")
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_simE_Pr
h_temp.SetFillColor(rt.kGreen-3)
h_temp.Draw("bar")
h_stack_Pr.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp.Add(h_trueE_Pr, h_simE_Pr, 1, -1)
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.7, 0.1, 0.9, 0.3)
legend.AddEntry(h_simE_Pr,"Good Reco","f")
legend.AddEntry(h_temp,"Improper Reco","f")
h_temp.GetListOfFunctions().Add(legend)
h_temp.Draw("bar")
h_stack_Pr.Add(h_temp)
legend.Draw()
h_stack_Pr.Draw("bar")
c_stack_Pr.Update()
c_stack_Pr.Write()

c_stack_Neu = rt.TCanvas("c_stack_Neu")
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_simE_Neu
h_temp.SetFillColor(rt.kGreen-3)
h_temp.Draw("bar")
h_stack_Neu.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp.Add(h_trueE_Neu, h_simE_Neu, 1, -1)
h_temp.SetFillColor(rt.kRed-7)
legend = rt.TLegend(0.7, 0.1, 0.9, 0.3)
legend.AddEntry(h_simE_Neu,"Good Reco","f")
legend.AddEntry(h_temp,"Improper Reco","f")
h_temp.GetListOfFunctions().Add(legend)
h_temp.Draw("bar")
h_stack_Neu.Add(h_temp)
legend.Draw()
h_stack_Neu.Draw("bar")
c_stack_Neu.Update()
c_stack_Neu.Write()

h_improper.Add(h_missing, h_misID, 1, 1)
h_improper.Add(h_secondary, 1)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_missing
h_temp.SetFillColor(rt.kRed-7)
h_temp.Draw("bar")
h_unconstructed_stack.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_misID
h_temp.SetFillColor(rt.kBlue-7)
h_temp.Draw("bar")
h_unconstructed_stack.Add(h_temp)
h_temp = rt.TH1F("h_temp","", 40,0,1200)
h_temp = h_secondary
h_temp.SetFillColor(rt.kYellow-7)
c_unconstructed_stack = rt.TCanvas("c_unconstructed_stack")
legend = rt.TLegend(0.7, 0.1, 0.9, 0.3)
legend.AddEntry(h_misID,"Wrong PDG","f")
legend.AddEntry(h_secondary,"Secondary Track","f")
legend.AddEntry(h_missing,"Missing Track","f")
h_temp_leg = h_temp.Clone("h_temp_leg")
h_temp_leg.GetListOfFunctions().Add(legend)
h_temp_leg.Draw("bar")
legend.Draw()
h_unconstructed_stack.Add(h_temp_leg)
h_unconstructed_stack.Draw("bar")
c_unconstructed_stack.Update()
#c_unconstructed_stack.Draw("bar")
c_unconstructed_stack.Write()


h_Eff_cc.Divide(h_simE_cc, h_trueE_cc)
h_Eff_Mu.Divide(h_simE_Mu, h_trueE_Mu)
h_Eff_Pr.Divide(h_simE_Pr, h_trueE_Pr)
h_Eff_Neu.Divide(h_simE_Neu, h_trueE_Neu)

#c_simE_cc = rt.TCanvas("c_simE_cc")
#h_simE_cc.Draw("bar")
#c_simE_cc.Update()
#c_simE_cc.Draw("bar")
#c_simE_cc.Write()
#
#c_trueE_cc = rt.TCanvas("c_trueE_cc")
#h_trueE_cc.Draw("bar")
#c_trueE_cc.Update()
#c_trueE_cc.Draw("bar")
#c_trueE_cc.Write()


c_Eff_cc = rt.TCanvas("c_Eff_cc")
h_Eff_cc.Draw()
c_Eff_cc.Update()
#c_Eff_cc.Draw("bar")
c_Eff_cc.Write()

c_recoComp_cc = rt.TCanvas("c_recoComp_cc")
h_recoComp_cc.Draw("bar")
c_recoComp_cc.Update()
#c_recoComp_cc.Draw("bar")
c_recoComp_cc.Write()

#c_simE_Mu = rt.TCanvas("c_simE_Mu")
#h_simE_Mu.Draw("bar")
#c_simE_Mu.Update()
#c_simE_Mu.Draw("bar")
#c_simE_Mu.Write()
#
#c_trueE_Mu = rt.TCanvas("c_trueE_Mu")
#h_trueE_Mu.Draw("bar")
#c_trueE_Mu.Update()
#c_trueE_Mu.Draw("bar")
#c_trueE_Mu.Write()

c_Eff_Mu = rt.TCanvas("c_Eff_Mu")
h_Eff_Mu.Draw()
c_Eff_Mu.Update()
#c_Eff_Mu.Draw("bar")
c_Eff_Mu.Write()

c_recoComp_Mu = rt.TCanvas("c_recoComp_Mu")
h_recoComp_Mu.Draw("bar")
c_recoComp_Mu.Update()
#c_recoComp_Mu.Draw("bar")
c_recoComp_Mu.Write()

#c_simE_Pr = rt.TCanvas("c_simE_Pr")
#h_simE_Pr.Draw("bar")
#c_simE_Pr.Update()
#c_simE_Pr.Draw("bar")
#c_simE_Pr.Write()
#
#c_trueE_Pr = rt.TCanvas("c_trueE_Pr")
#h_trueE_Pr.Draw("bar")
#c_trueE_Pr.Update()
#c_trueE_Pr.Draw("bar")
#c_trueE_Pr.Write()

c_Eff_Pr = rt.TCanvas("c_Eff_Pr")
h_Eff_Pr.Draw()
c_Eff_Pr.Update()
#c_Eff_Pr.Draw("bar")
c_Eff_Pr.Write()

c_recoComp_Pr = rt.TCanvas("c_recoComp_Pr")
h_recoComp_Pr.Draw("bar")
c_recoComp_Pr.Update()
#c_recoComp_Pr.Draw("bar")
c_recoComp_Pr.Write()

#c_simE_Neu = rt.TCanvas("c_simE_Neu")
#h_simE_Neu.Draw("bar")
#c_simE_Neu.Update()
#c_simE_Neu.Draw("bar")
#c_simE_Neu.Write()
#
#c_trueE_Neu = rt.TCanvas("c_trueE_Neu")
#h_trueE_Neu.Draw("bar")
#c_trueE_Neu.Update()
#c_trueE_Neu.Draw("bar")
#c_trueE_Neu.Write()

c_Eff_Neu = rt.TCanvas("c_Eff_Neu")
h_Eff_Neu.Draw()
c_Eff_Neu.Update()
#c_Eff_Neu.Draw("bar")
c_Eff_Neu.Write()

c_recoComp_Neu = rt.TCanvas("c_recoComp_Neu")
h_recoComp_Neu.Draw("bar")
c_recoComp_Neu.Update()
#c_recoComp_Neu.Draw("bar")
c_recoComp_Neu.Write()

c_missing = rt.TCanvas("c_missing")
h_missing.Draw("bar")
c_missing.Update()
#c_missing.Draw("bar")
c_missing.Write()

c_misID = rt.TCanvas("c_misID")
h_misID.Draw("bar")
c_misID.Update()
#c_misID.Draw("bar")
c_misID.Write()

c_secondary = rt.TCanvas("c_secondary")
h_secondary.Draw("bar")
c_secondary.Update()
c_secondary.Draw("bar")
c_secondary.Write()

c_improper = rt.TCanvas("c_improper")
h_improper.Draw("bar")
c_improper.Update()
#c_improper.Draw("bar")
c_improper.Write()
