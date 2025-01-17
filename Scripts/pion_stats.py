# pion_stats.py
# Author: Andy Zhang
# Purpose: Categorize pions by truth info, mainly trajectory history, colinearity with other primaries, and secondary particles.
# Data Created: July 2nd, 2024
# Last Modified: July 31st, 2024


import ROOT as rt
import numpy as np
from ROOT import gROOT, addressof
from array import array
# f1 = rt.TFile("inputFile.root")
f1 = rt.TFile("make_dlgen2_flat_ntuples_reco_v2me06_gen2ntuple_v5_output.root")
f_out = rt.TFile("pion_stats.root", "RECREATE")
pion_tree = rt.TTree("pion_tree", "pion_tree")
potTree = rt.TTree("potTree", "potTree")
t1 = f1.Get("EventTree")
p1 = f1.Get("potTree")

# Recursive function to find ancestor particle
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

def mom_angle(mv1, mv2):
    magni1 = np.linalg.norm(mv1)
    magni2 = np.linalg.norm(mv2)
    angle = (np.dot(mv1, mv2) / (magni1*magni2))
    return angle
        
class pion_container:
    def __init__(self, true_trackID, true_pionE, true_endE, xsecWeight, Contained, EventID, CCNC, init_pos, end_pos):
        self.true_trackID = true_trackID
        self.true_pionE = true_pionE
        self.true_endE = true_endE
        self.xsecWeight = xsecWeight
        self.contained = Contained
        self.EventID = EventID
        self.CCNC = CCNC
        self.init_pos = init_pos
        self.end_pos = end_pos
        self.colinear_ang = None
        self.colinear_TID = -1
        self.secondary_TIDs = []
        self.secondary_PDGs = []
        self.secondary_Source = []
        self.secondary_init_pos = []
        self.nsecondaries = 0
        self.reco_status = 0
        self.track_length = np.abs(np.linalg.norm(np.array(init_pos) - np.array(end_pos)))
    def set_colinear(self, colinear_ang, colinear_TID):
        self.colinear_ang = colinear_ang
        self.colinear_TID = colinear_TID
    def add_secondary(self, TID, PDG, Source, init_pos):
        self.secondary_TIDs.append(TID)
        self.secondary_PDGs.append(PDG)
        self.secondary_Source.append(Source)
        self.secondary_init_pos.append(init_pos)
        self.nsecondaries += 1
    def set_reco_status(self, reco_status):
        self.reco_status = reco_status
    
        

# Set up tree branches
p_TID = array('i', [0])
pion_tree.Branch('p_TID', p_TID, 'p_TID/I')
p_Tracklength = array('f', [0])
pion_tree.Branch('p_Tracklength', p_Tracklength, 'p_Tracklength/F')
p_TrueE = array('f', [0,0,0,0])
pion_tree.Branch('p_TrueE', p_TrueE, 'p_TrueE[4]/F')
p_Contained = array('i', [0])
pion_tree.Branch('p_Contained', p_Contained, 'p_Contained/I')
nSecondaries = array('i', [0])
pion_tree.Branch('nSecondaries', nSecondaries, 'nSecondaries/I')

p_EventID = array('i', [0])
pion_tree.Branch('p_EventID', p_EventID, 'p_EventID/I')
p_Run = array('i', [0])
pion_tree.Branch('p_Run', p_Run, 'p_Run/I')
p_Subrun = array('i', [0])
pion_tree.Branch('p_Subrun', p_Subrun, 'p_Subrun/I')
p_FileID = array('i', [0])
pion_tree.Branch('p_FileID', p_FileID, 'p_FileID/I')
p_EventWeight = array('f', [0])
pion_tree.Branch('p_EventWeight', p_EventWeight, 'p_EventWeight/F')

p_CCNC = array('i', [0])
pion_tree.Branch('p_CCNC', p_CCNC, 'p_CCNC/I')
p_colinearTID = array('i', [0])
pion_tree.Branch('p_colinearTID', p_colinearTID, 'p_colinearTID/I')
p_colinearAng = array('f', [0])
pion_tree.Branch('p_colinearAng', p_colinearAng, 'p_colinearTID/F')
p_History = array('i', [0])
pion_tree.Branch('p_History', p_History, 'p_History/I')
p_Dist_to_interaction = array('f', [0])
pion_tree.Branch('p_Dist_to_interaction', p_Dist_to_interaction, 'p_Dist_to_interaction/F')
p_Interactions = array('i', [0])
pion_tree.Branch('p_Interactions', p_Interactions, 'p_Interactions/I')

max_sec_parts = 50
sec_TruePDG = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TruePDG', sec_TruePDG, 'sec_TruePDG[nSecondaries]/I')
sec_TrueTID = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TrueTID', sec_TrueTID, 'sec_TrueTID[nSecondaries]/I')
sec_TrueSource = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TrueSource', sec_TrueSource, 'sec_TrueSource[nSecondaries]/I')

p_Reco_status = array('i', [0])
pion_tree.Branch('p_Reco_status', p_Reco_status, 'p_Reco_status/I')

totPOT = array('f', [0])
potTree.Branch('totPOT', totPOT, 'totPOT/F')
totGoodPOT = array('f', [0])
potTree.Branch('totGoodPOT', totGoodPOT, 'totGoodPOT/F')

for entry in p1:
    POT = entry.totPOT
    GoodPOT = entry.totGoodPOT

    totPOT[0] = POT
    totGoodPOT[0] = GoodPOT
    potTree.Fill()

prim_sec_counter = [0,0]
# Loop over all entries
for entry in t1:
    xsecWeight = entry.xsecWeight
    trueNuPDG = entry.trueNuPDG
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
    trueSimPartContained = entry.trueSimPartContained
    trueSimPartX = entry.trueSimPartX
    trueSimPartY = entry.trueSimPartY
    trueSimPartZ = entry.trueSimPartZ
    trueSimPartEndX = entry.trueSimPartEndX
    trueSimPartEndY = entry.trueSimPartEndY
    trueSimPartEndZ = entry.trueSimPartEndZ
    event = entry.event
    run = entry.run
    subrun = entry.subrun
    fileid = entry.fileid
    trueSimPartEndE = entry.trueSimPartEndE
    trueSimPartEndPx = entry.trueSimPartEndPx
    trueSimPartEndPy = entry.trueSimPartEndPy
    trueSimPartEndPz = entry.trueSimPartEndPz

    nTracks = entry.nTracks
    trackTrueTID = entry.trackTrueTID
    trackPID = entry.trackPID
    trackIsSecondary = entry.trackIsSecondary

    trueVtxX = entry.trueVtxX
    trueVtxY = entry.trueVtxY
    trueVtxZ = entry.trueVtxZ
    foundVertex = entry.foundVertex
    vtxX = entry.vtxX
    vtxY = entry.vtxY
    vtxZ = entry.vtxZ

    # Skip the entry if the vertex wasn't reconstructed or was too far away
    trueVtx = np.array([trueVtxX, trueVtxY, trueVtxZ])
    recoVtx = np.array([vtxX, vtxY, vtxZ])
    vtx_dist = np.linalg.norm(trueVtx - recoVtx) 
    if (foundVertex == 0 or vtx_dist > 3):
        continue
    
    # Find pion truth Data
    pion_mv = 0
    event_pions = []
    ancestor_index = []
    for i in range(0, nTrueSimParts):
        if trueSimPartProcess[i] == 0:
            if (np.abs(trueSimPartPDG[i] == 211)):
                init_pos = [trueSimPartX[i], trueSimPartY[i], trueSimPartZ[i]]
                end_pos = [trueSimPartEndX[i], trueSimPartEndY[i], trueSimPartEndZ[i]]
                TrueE = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i], trueSimPartE[i]]
                EndE = [trueSimPartEndPx[i], trueSimPartEndPy[i], trueSimPartEndPz[i], trueSimPartEndE[i]]
                Contained = trueSimPartContained[i]
                TID = trueSimPartTID[i]
                EventID = [fileid,run,subrun,event]
                EventWeight = xsecWeight
                CCNC = trueNuCCNC
                pion_mv = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i]]
                pion = pion_container(TID, TrueE, EndE, EventWeight, Contained, EventID, CCNC, init_pos, end_pos)
                event_pions.append(pion)
                # print("pion found")
        
        for k in range(0,nTracks):
            if (trueSimPartTID[i] == trackTrueTID[k]):
                if(trueSimPartProcess[i] == 0):
                    ancestor_index.append([prim_ancestor(trueSimPartTID, trackTrueTID[k]),0])
                else:
                    ancestor_index.append([prim_ancestor(trueSimPartTID, trackTrueTID[k]),1])

    # Find pion-secondaries and primaries truth data
    most_colinear = [0,0]
    for j in range(0,len(event_pions)):
        for i in range(0, nTrueSimParts):
            # Add a secondary particle to the pion container
            if (event_pions[j].true_trackID == trueSimPartMID[i] and trueSimPartProcess[i] != 0):
                if (np.abs(trueSimPartPDG[i]) != 211):
                    init_pos = [trueSimPartX[i], trueSimPartY[i], trueSimPartZ[i]]
                    event_pions[j].add_secondary(trueSimPartTID[i], trueSimPartPDG[i], trueSimPartProcess[i], init_pos)
            
            # Calculate colinear angle to the pion of current primary
            
            ##### FIX COLINEARITY CALCULATION TO USE THE TRACK SEGMENT BEFORE ANY INTERACTION
            if (trueSimPartProcess[i] == 0):
                if (trueSimPartTID[i] != event_pions[j].true_trackID):
                    pion_mv = [event_pions[j].true_pionE[0], event_pions[j].true_pionE[1], event_pions[j].true_pionE[2]]
                    new_mv = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i]]
                    colinearity = mom_angle(pion_mv, new_mv)
                    if (event_pions[j].colinear_ang == None):
                        event_pions[j].set_colinear(colinearity, trueSimPartTID[i])
                    elif (np.abs(colinearity) >= np.abs(event_pions[j].colinear_ang)):
                        event_pions[j].set_colinear(colinearity, trueSimPartTID[i])

        TID_found = False
        reco_found = False
        sec_index = []
        prim_index = []
        for k in range(0, nTracks):
            if(trackIsSecondary[k] != 1):
                if(event_pions[j].true_trackID == trackTrueTID[k]):
                    TID_found = True
                    if(np.abs(trackPID[k] == 211)):
                        reco_found = True

                track_ancestor = prim_ancestor(trueSimPartTID, trackTrueTID[k])
                prim_index.append(track_ancestor)
            else:
                track_ancestor = prim_ancestor(trueSimPartTID, trackTrueTID[k])
                sec_index.append(track_ancestor)
                
            # Fill histograms according to pion recostruction type
        if(TID_found == True):
            if(reco_found == True):
                event_pions[j].set_reco_status(0)
            else: # Track was found, but incorrect PDG
                event_pions[j].set_reco_status(1)
        else:
            sec_found = False
            prim_found = False
            for i in range(0,len(ancestor_index)):
                if (len(ancestor_index[i]) == 0 or ancestor_index[i][0] == 999):
                    continue
                if (ancestor_index[i][1] == 0 and trueSimPartMID[ancestor_index[i][0]] == event_pions[j].true_trackID):
                    prim_found = True
                if (ancestor_index[i][1] == 1 and trueSimPartTID[ancestor_index[i][0]] == event_pions[j].true_trackID):
                    sec_found = True

            # for i in range(0, len(sec_index)):
            #     if(sec_index[i] == 999):
            #         continue
            #     # Track not found, but a secondary of the pion was found
            #     if(trueSimPartTID[sec_index[i]] == event_pions[j].true_trackID):
            #         sec_found = True
            #         prim_sec_counter[0] += 1
            #         break
            # for i in range(0, len(prim_index)):
            #     if(prim_index[i] == 999):
            #         continue
            #     # Track not found, but a secondary of the pion was found
            #     if(trueSimPartMID[prim_index[i]] == event_pions[j].true_trackID):
            #         prim_found = True
            #         prim_sec_counter[1] += 1
            #         break
            # Track not found, secondaries found
            if(sec_found == True):
                event_pions[j].set_reco_status(2)
                prim_sec_counter[1] += 1
            elif (prim_found == True): # Track not found, secondaries not found
                event_pions[j].set_reco_status(3)
                prim_sec_counter[0] += 1
            else: 
                event_pions[j].set_reco_status(4)

    # Fill branch info
    for i in range(0,len(event_pions)):
        p_TID[0] = event_pions[i].true_trackID
        p_Tracklength = event_pions[i].track_length
        p_TrueE[0] = event_pions[i].true_pionE[0]
        p_TrueE[1] = event_pions[i].true_pionE[1]
        p_TrueE[2] = event_pions[i].true_pionE[2]
        p_TrueE[3] = event_pions[i].true_pionE[3]
        p_Contained[0] = event_pions[i].contained

        p_FileID[0] = event_pions[i].EventID[0]
        p_Run[0] = event_pions[i].EventID[1]
        p_Subrun[0] = event_pions[i].EventID[2]
        p_EventID[0] = event_pions[i].EventID[3]
        p_EventWeight[0] = event_pions[i].xsecWeight

        p_CCNC[0] = event_pions[i].CCNC
        p_colinearTID[0] = event_pions[i].colinear_TID
        colinearity = event_pions[i].colinear_ang
        if colinearity == None:   
            p_colinearAng[0] = -1.1
        else:
            p_colinearAng[0] = colinearity

        p_Reco_status[0] = event_pions[i].reco_status
        
        
        nSecondaries[0] = event_pions[i].nsecondaries
        secondary_mode = [0,0]

        # Set up secondary branches dependent on nSecondaries
        for j in range(0,nSecondaries[0]):
            sec_TruePDG[j] = event_pions[i].secondary_PDGs[j]
            sec_TrueTID[j] = event_pions[i].secondary_TIDs[j]
            sec_TrueSource[j] = event_pions[i].secondary_Source[j]
            if (event_pions[i].secondary_Source[j] == 1):
                secondary_mode[0] = 1 # decay
            elif (event_pions[i].secondary_Source[j] == 2):
                secondary_mode[1] = 1 # other processes

        # record what happened to the pion; 0 for un-contained, 1 for contained and ranged out, 2 for decay in flight, 3 for secondary interactions
        p_History[0] = 4
        p_Interactions[0] = 0
        KE = event_pions[i].true_endE[3] - np.sqrt(event_pions[i].true_endE[3]**2 
        - (event_pions[i].true_endE[0]**2 + event_pions[i].true_endE[1]**2 + event_pions[i].true_endE[2]**2))
        if (secondary_mode[1] == 1):
            p_History[0] = 2
            # print("SEC INTERACT")
        elif (secondary_mode[0] == 1 and KE > 0.010):
            p_History[0] = 3
            # print("DECAY IN FLIGHT")
        elif(KE <= 0.010): # 10 keV
            p_History[0] = 1
            # print("RANGE OUT")
        elif (event_pions[i].contained == 0):
            p_History[0] = 0
            # print("UNCONTAINED")

        if (secondary_mode[1] == 1):
            p_Interactions[0] = 1
            # print("SEC INTERACT")
        if (secondary_mode != 0):
            int_dist = 10000
            for j in range(0,nSecondaries[0]):
                test_length = np.linalg.norm(np.array(event_pions[i].init_pos) 
                - np.array(event_pions[i].secondary_init_pos[j]))
                if test_length < int_dist:
                    int_dist = test_length
            if (event_pions[i].track_length < int_dist):
                int_dist = event_pions[i].track_length
            p_Dist_to_interaction[0] = int_dist
        if (p_History[0] == 4):
            print("ERROR: EXCEPTION FOR PION MODE, KE is ", KE)
        pion_tree.Fill()

f_out.cd()
pion_tree.Write("",rt.TObject.kOverwrite)
potTree.Write("",rt.TObject.kOverwrite)
f_out.Close()

print(prim_sec_counter)
# Also add extra information on the pion secondaries' histories, i.e. if they have secondary interactions and decays. Take note of
