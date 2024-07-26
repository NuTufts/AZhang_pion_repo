import ROOT as rt
import numpy as np
from ROOT import gROOT, AddressOf
f1 = rt.TFile("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_nu_overlay.root")
f2 = rt.TFile("dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_intrinsic_nue_overlay.root")
t1 = f1.Get("EventTree")
t2 = f2.Get("EventTree")
p1 = f1.Get("potTree")
p2 = f2.Get("potTree")

class pion_container:
    def __init__(self, true_trackID, top_recoComp, true_pionE, xsecWeight):
        self.true_trackID = true_trackID
        self.reco_found = False
        self.top_recoComp = top_recoComp
        self.true_pionE = true_pionE
        self.xsecWeight = xsecWeight
    def toggle_reco_found(self):
        self.reco_found = True
    def set_top_recoComp(self, trackTrueComp):
        self.top_recoComp = trackTrueComp

t_pot = 0.
totGoodPOT_1 = 0.
totGoodPOT_2 = 0.

for entry in p1:
    t_pot += entry.totGoodPOT
for entry in p2:
    t_pot += entry.totGoodPOT

h_simE_cc = rt.TH1F("h_simE_cc","Number of Simulated Pions by Energy (CC);True Pion Energy (MeV);# of Pions per 6.67e20 POT", 40,0,4000)

h_trueE_cc = rt.TH1F("h_trueE_cc","Number of Truth-Matched Pions by Energy (CC);True Pion Energy (MeV);# of Pions per 6.67e20 POT", 40,0,4000)

h_Eff_cc = rt.TH1F("h_Eff_cc","Pion Reconstruction Efficiency by Energy (CC);True Pion Energy (MeV);Reconstruction Efficiency", 40,0,4000)

h_simE_nc = rt.TH1F("h_simE_nc","Number of Simulated Pions by Energy (NC);True Pion Energy (MeV);# of Pions per 6.67e20 POT", 40,0,4000)

h_trueE_nc = rt.TH1F("h_trueE_nc","Number of Truth-Matched Pions by Energy (NC);True Pion Energy (MeV);# of Pions per 6.67e20 POT", 40,0,4000)

h_Eff_nc = rt.TH1F("h_Eff_nc","Pion Reconstruction Efficiency by Energy (NC);True Pion Energy (MeV);Reconstruction Efficiency", 40,0,4000)

h_recoComp_cc = rt.TH1F("h_recoComp_cc","Greatest Track Completeness of Truth-Matched Pion Events (CC);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

h_recoComp_nc = rt.TH1F("h_recoComp_nc","Greatest Track Completeness of Truth-Matched Pion Events (NC);Reconstructed Track Completeness (%);# of Truth-Matched Pion Events", 20,0,1)

cc = []
nc = []
for entry in t1:
    xsecWeight = entry.xsecWeight
    trueNuPDG = entry.trueNuPDG
    nTracks = entry.nTracks
    trackTrueTID = entry.trackTrueTID
    nTrueSimParts = entry.nTrueSimParts
    trueSimPartTID = entry.trueSimPartTID
    trueSimPartE = entry.trueSimPartE
    trueSimPartPDG = entry.trueSimPartPDG
    trueSimPartProcess = entry.trueSimPartProcess
    trueNuCCNC = entry.trueNuCCNC
    trackPID = entry.trackPID
    trackTrueComp = entry.trackTrueComp
    
   #if(trueNuPDG == 12 or trueNuPDG == -12):
    if(trueNuCCNC == 1):
        for j in range(0, nTrueSimParts):
            if(trueSimPartProcess[j] == 0):
                if(trueSimPartPDG[j] == 211 or trueSimPartPDG[j] == -211):
                    h_trueE_cc.Fill(trueSimPartE[j], xsecWeight)
                    temp = pion_container(trueSimPartTID[j], 0, trueSimPartE[j], xsecWeight)
                    cc.append(temp)

        for j in range(0, len(cc)):
            for k in range(0, nTracks):
                if(cc[j].true_trackID == trackTrueTID[k]):
                    if(trackPID[k] == 211 or trackPID[k] == -211):
                        cc[j].toggle_reco_found()
                        if(trackTrueComp[k] > cc[j].top_recoComp):
                            cc[j].set_top_recoComp(trackTrueComp[k])
                        h_simE_cc.Fill(cc[j].true_pionE, cc[j].xsecWeight)
                        break
    else:
        for j in range(0, nTrueSimParts):
            if(trueSimPartProcess[j] == 0):
                if(trueSimPartPDG[j] == 211 or trueSimPartPDG[j] == -211):
                    h_trueE_nc.Fill(trueSimPartE[j], xsecWeight)
                    temp = pion_container(trueSimPartTID[j], 0, trueSimPartE[j], xsecWeight)
                    nc.append(temp)

        for j in range(0, len(nc)):
            for k in range(0, nTracks):
                if(nc[j].true_trackID == trackTrueTID[k]):
                    if(trackPID[k] == 211 or trackPID[k] == -211):
                        nc[j].toggle_reco_found()
                        if(trackTrueComp[k] > nc[j].top_recoComp):
                            nc[j].set_top_recoComp(trackTrueComp[k])
                        h_simE_nc.Fill(nc[j].true_pionE, nc[j].xsecWeight)
                        break
                        
for entry in t1:
    xsecWeight = entry.xsecWeight
    trueNuPDG = entry.trueNuPDG
    nTracks = entry.nTracks
    trackTrueTID = entry.trackTrueTID
    nTrueSimParts = entry.nTrueSimParts
    trueSimPartTID = entry.trueSimPartTID
    trueSimPartE = entry.trueSimPartE
    trueSimPartPDG = entry.trueSimPartPDG
    trueSimPartProcess = entry.trueSimPartProcess
    trueNuCCNC = entry.trueNuCCNC
    trackPID = entry.trackPID
    trackTrueComp = entry.trackTrueComp
    
   #if(trueNuPDG == 12 or trueNuPDG == -12):
    if(trueNuCCNC == 1):
        for j in range(0, nTrueSimParts):
            if(trueSimPartProcess[j] == 0):
                if(trueSimPartPDG[j] == 211 or trueSimPartPDG[j] == -211):
                    h_trueE_cc.Fill(trueSimPartE[j], xsecWeight)
                    temp = pion_container(trueSimPartTID[j], 0, trueSimPartE[j], xsecWeight)
                    cc.append(temp)

        for j in range(0, len(cc)):
            for k in range(0, nTracks):
                if(cc[j].true_trackID == trackTrueTID[k]):
                    if(trackPID[k] == 211 or trackPID[k] == -211):
                        cc[j].toggle_reco_found()
                        if(trackTrueComp[k] > cc[j].top_recoComp):
                            cc[j].set_top_recoComp(trackTrueComp[k])
                        h_simE_cc.Fill(cc[j].true_pionE, cc[j].xsecWeight)
                        break
    else:
        for j in range(0, nTrueSimParts):
            if(trueSimPartProcess[j] == 0):
                if(trueSimPartPDG[j] == 211 or trueSimPartPDG[j] == -211):
                    h_trueE_nc.Fill(trueSimPartE[j], xsecWeight)
                    temp = pion_container(trueSimPartTID[j], 0, trueSimPartE[j], xsecWeight)
                    nc.append(temp)

        for j in range(0, len(nc)):
            for k in range(0, nTracks):
                if(nc[j].true_trackID == trackTrueTID[k]):
                    if(trackPID[k] == 211 or trackPID[k] == -211):
                        nc[j].toggle_reco_found()
                        if(trackTrueComp[k] > nc[j].top_recoComp):
                            nc[j].set_top_recoComp(trackTrueComp[k])
                        h_simE_nc.Fill(nc[j].true_pionE, nc[j].xsecWeight)
                        break
                        

for i in range(0, len(cc)):
    if(cc[i].reco_found == True):
        h_recoComp_cc.Fill(cc[i].top_recoComp)
for i in range(0, len(nc)):
    if(nc[i].reco_found == True):
        h_recoComp_cc.Fill(nc[i].top_recoComp)
        
h_simE_cc.Scale(6.67e20/t_pot)
h_trueE_cc.Scale(6.67e20/t_pot)
h_simE_nc.Scale(6.67e20/t_pot)
h_trueE_nc.Scale(6.67e20/t_pot)
h_Eff_cc.Divide(h_trueE_cc, h_simE_cc)
h_Eff_nc.Divide(h_trueE_nc, h_simE_nc)
h_simE_cc.Draw("EHIST")
h_trueE_cc.Draw("EHIST")
h_Eff_cc.Draw("EHIST")
h_recoComp_cc.Draw("EHIST")
h_simE_nc.Draw("EHIST")
h_trueE_nc.Draw("EHIST")
h_Eff_nc.Draw("EHIST")
h_recoComp_nc.Draw("EHIST")

f_out = rt.TFile("Pi_Efficiency.root","RECREATE")
h_simE_cc.Write()
h_trueE_cc.Write()
h_Eff_cc.Write()
h_recoComp_cc.Write()
h_simE_nc.Write()
h_trueE_nc.Write()
h_Eff_nc.Write()
h_recoComp_nc.Write()
