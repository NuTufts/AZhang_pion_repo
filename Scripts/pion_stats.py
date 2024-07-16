# pion_stats.py
# Author: Andy Zhang
# Purpose: Categorize pions by truth info, mainly trajectory history, colinearity with other primaries, and secondary particles.
# Data Created: July 2nd, 2024
# Last Modified: July 11th, 2024


import ROOT as rt
import numpy as np
from ROOT import gROOT, addressof
from array import array
f1 = rt.TFile("inputFile.root")
f_out = rt.TFile("pion_stats.root", "RECREATE")
pion_tree = rt.TTree("pion_tree", "pion_tree")
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
    angle = np.absolute(np.arccos(((np.dot(mv1, mv2) / magni1*magni2)) % 1))
    return angle
        
class pion_container:
    def __init__(self, true_trackID, true_pionE, true_endE, xsecWeight, Contained, EventID, CCNC, init_pos):
        self.true_trackID = true_trackID
        self.true_pionE = true_pionE
        self.true_endE = true_endE
        self.xsecWeight = xsecWeight
        self.contained = Contained
        self.EventID = EventID
        self.CCNC = CCNC
        self.init_pos = init_pos
        self.colinear_ang = 0
        self.colinear_TID = 0
        self.secondary_TIDs = []
        self.secondary_PDGs = []
        self.secondary_Source = []
        self.secondary_init_pos = []
        self.nsecondaries = 0
        self.tracklength = 0
    def set_colinear(self, colinear_ang, colinear_TID):
        self.colinear_ang = colinear_ang
        self.colinear_TID = colinear_TID
    def add_secondary(self, TID, PDG, Source, init_pos):
        self.secondary_TIDs.append(TID)
        self.secondary_PDGs.append(PDG)
        self.secondary_Source.append(Source)
        self.secondary_init_pos.append(init_pos)
        self.nsecondaries += 1
    def set_tracklength(self, tracklength):
        self.tracklength = tracklength
    
        

# Set up tree branches
p_TID = array('i', [0])
pion_tree.Branch('p_TID', p_TID, 'p_TID/I')
p_Tracklength = array('f', [0])
pion_tree.Branch('p_Tracklength', p_Tracklength, 'p_Tracklength/F')
p_TrueE = array('f', [0,0,0,0])
pion_tree.Branch('p_TrueE', p_TrueE, 'p_TrueE/F')
p_Contained = array('i', [0])
pion_tree.Branch('p_Contained', p_Contained, 'p_Contained/I')
nSecondaries = array('i', [0])
pion_tree.Branch('nSecondaries', nSecondaries, 'nSecondaries/I')
p_EventID = array('i', [0])
pion_tree.Branch('p_EventID', p_EventID, 'p_EventID/I')
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
p_Dist_to_interaction = array('i', [0])
pion_tree.Branch('p_Dist_to_interaction', p_Dist_to_interaction, 'p_Dist_to_interaction/I')
max_sec_parts = 50
sec_TruePDG = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TruePDG', sec_TruePDG, 'sec_TruePDG[nSecondaries]/I')
sec_TrueTID = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TrueTID', sec_TrueTID, 'sec_TrueTID[nSecondaries]/I')
sec_TrueSource = array('i', [0]*max_sec_parts)
pion_tree.Branch('sec_TrueSource', sec_TrueSource, 'sec_TrueSource[nSecondaries]/I')


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
    event = entry.event
    # trueSimPartEndE = entry.trueSimPartEndE
    # trueSimPartEndPx = entry.trueSimPartEndPx
    # trueSimPartEndPy = entry.trueSimPartEndPy
    # trueSimPartEndPz = entry.trueSimPartEndPz

    # Find pion truth Data
    pion_mv = 0
    event_pions = []
    for i in range(0, nTrueSimParts):
        if trueSimPartProcess[i] == 0:
            if (trueSimPartPDG[i] == 211 or trueSimPartPDG[i] == -211):
                init_pos = [trueSimPartX, trueSimPartY, trueSimPartZ]
                TrueE = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i], trueSimPartE[i]]
                # EndE = [trueSimPartEndPx[i], trueSimPartEndPy[i], trueSimPartEndPz[i], trueSimPartEndE[i]]
                Contained = trueSimPartContained[i]
                TID = trueSimPartTID[i]
                EventID = event
                EventWeight = xsecWeight
                CCNC = trueNuCCNC
                pion_mv = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i]]
                pion = pion_container(TID, TrueE, EventWeight, Contained, EventID, CCNC, init_pos)
                # pion = pion_container(TID, TrueE, EndE, EventWeight, Contained, EventID, CCNC, init_pos)
                event_pions.append(pion)
#                print("pion found")

    # Find pion-secondaries and primaries truth data
    most_colinear = [0,0]
    for j in range(0,len(event_pions)):
        for i in range(0, nTrueSimParts):
            index = prim_ancestor(trueSimPartTID, trueSimPartTID[i])
            if index == 999:
                continue
            # Add a secondary particle to the pion container
            if (event_pions[j].true_trackID == trueSimPartMID[i] and trueSimPartProcess[i] != 0):
                init_pos = [trueSimPartX, trueSimPartY, trueSimPartZ]
                event_pions[j].add_secondary(trueSimPartTID[i], trueSimPartPDG[i], trueSimPartProcess[i], init_pos)
            
            # Calculate colinear angle to the pion of current primary
            
            ##### FIX COLINEARITY CALCULATION TO USE THE TRACK SEGMENT BEFORE ANY INTERACTION
            if (trueSimPartProcess[i] == 0):
                new_mv = [trueSimPartPx[i], trueSimPartPy[i], trueSimPartPz[i]]
                colinearity = mom_angle(pion_mv, new_mv)
                if (colinearity <= event_pions[j].colinear_ang):
                    event_pions[j].set_colinear(colinearity, trueSimPartTID[i])
    
    # Fill branch info
    for i in range(0,len(event_pions)):
        p_TID[0] = event_pions[i].true_trackID
        p_Tracklength[0] = event_pions[i].tracklength
        p_TrueE[0] = event_pions[i].true_pionE[0]
        p_TrueE[1] = event_pions[i].true_pionE[1]
        p_TrueE[2] = event_pions[i].true_pionE[2]
        p_TrueE[3] = event_pions[i].true_pionE[3]
        p_Contained[0] = event_pions[i].contained
        p_EventID[0] = event_pions[i].EventID
        p_EventWeight[0] = event_pions[i].xsecWeight
        p_CCNC[0] = event_pions[i].CCNC
        p_colinearTID[0] = event_pions[i].colinear_TID
        p_colinearAng[0] = event_pions[i].colinear_ang
        
        nSecondaries[0] = event_pions[i].nsecondaries
        secondary_mode = 0
#        print("pion TID: " + str(p_TID[0]))

        # Set up secondary branches dependent on nSecondaries
        for j in range(0,nSecondaries[0]):
            sec_TruePDG[j] = event_pions[i].secondary_PDGs[j]
            sec_TrueTID[j] = event_pions[i].secondary_TIDs[j]
            sec_TrueSource[j] = event_pions[i].secondary_Source[j]
            if (event_pions[i].secondary_Source[j] == 1):
                secondary_mode = 1
            elif (event_pions[i].secondary_Source[j] == 2):
                secondary_mode = 2

        # record what happened to the pion; 0 for un-contained, 1 for contained and ranged out, 2 for decay in flight, 3 for secondary interactions
        p_History[0] = 999
        if (event_pions[i].contained == 0):
            p_History[0] = 0
            print("UNCONTAINED")
        elif (np.linalg.norm(np.array(event_pions[i].true_EndE)) == 0):
            p_History[0] = 1
            print("RANGE OUT")
        elif (secondary_mode == 1):
            p_History[0] = 2
            tracklength = np.linalg.norm(np.array(event_pions[i].init_pos) - np.array(event_pions[i].secondary_init_pos[0]))
            p_Tracklength[0] = tracklength
            print("DECAY")
        elif (secondary_mode == 2):
            p_History[0] = 3
            tracklength = np.linalg.norm(np.array(event_pions[i].init_pos) - np.array(event_pions[i].secondary_init_pos[0]))
            p_Tracklength[0] = tracklength
            print("SEC INTERACT")
        if (p_History[0] == 999):
            print("ERROR: EXCEPTION FOR PION MODE")
        pion_tree.Fill()

f_out.cd()
pion_tree.Write("",rt.TObject.kOverwrite)
f_out.Close()

# Also add extra information on the pion secondaries' histories, i.e. if they have secondary interactions and decays. Take note of
