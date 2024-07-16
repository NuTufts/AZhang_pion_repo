
import os,sys,argparse

import ROOT as rt

from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow

from math import sqrt as sqrt
from math import acos as acos
from math import pi
from math import isinf

from array import array
import numpy as np

from event_weighting.event_weight_helper import SumPOT, Weights
from helpers.larflowreco_ana_funcs import *
from helpers.pionEnergyEstimator import pionRange2T

import torch
from torch import nn
from torch.utils.data import DataLoader
import torchvision.transforms as transforms

parser = argparse.ArgumentParser("Make Flat NTuples for DLGen2 Analyses")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list or merged_dlreco file for single input")
parser.add_argument("-w", "--weightfile", type=str, default="none", help="weights file (pickled python dict)")
parser.add_argument("-m", "--model_path", type=str, required=True, help="path to prong CNN checkpoint file")
parser.add_argument("-d", "--device", type=str, default="cpu", help="gpu/cpu device")
parser.add_argument("-mc", "--isMC", help="running over MC input", action="store_true")
parser.add_argument("-ana", "--dlana_input", help="using merged_dlana input files", action="store_true")
parser.add_argument("-o", "--outfile", type=str, default="dlgen2_flat_ntuple.root", help="output file name")
parser.add_argument("-nkp","--noKeypoints", action="store_true", help="don't save keypoint info")
parser.add_argument("--ignoreWeights", action="store_true", help="don't look up xsec weights, set to 1 and process all MC events (default: lookup xsec weights and exit with error if not found)")
parser.add_argument("--skipNoWeightEvts", action="store_true", help="skip MC events if we can't find xsec weights but continue processing (default: exit with error)")
parser.add_argument("--multiGPU", action="store_true", help="use multiple GPUs")
args = parser.parse_args()

sys.path.append(args.model_path[:args.model_path.find("/checkpoints")])
from models_instanceNorm_reco_2chan_quadTask import ResBlock, ResNet34
from normalization_constants import mean, std

if args.isMC and args.weightfile=="none" and not args.ignoreWeights:
  sys.exit("Must supply weight file for MC input. Exiting...")

reco2Tag = "merged_dlreco_"
if args.dlana_input:
  reco2Tag = "merged_dlana_"

if ".root" in args.truth and len(args.files) == 1 and ".root" in args.files[0]:
  files = [ [args.files[0], args.truth] ]
else:
  if len(args.files) == 1 and ".txt" in args.files[0]:
    larflowfiles = []
    with open(args.files[0], "r") as larflowlist:
      for line in larflowlist:
        larflowfiles.append(line.replace("\n",""))
    files = getFiles(reco2Tag, larflowfiles, args.truth)
  else:
    files = getFiles(reco2Tag, args.files, args.truth)


def addClusterCharge(iolcv, cluster, vertexPixels, vertexCharge, threshold):
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  clusterPixels = []
  clusterCharge = 0.
  for hit in cluster:
    for p in range(3):
      row = (hit.tick - 2400)//6
      pixel = [ p, hit.tick, hit.targetwire[p] ]
      pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
      if pixVal < threshold:
        continue
      if pixel not in clusterPixels:
        clusterPixels.append(pixel)
        clusterCharge += pixVal
      if pixel not in vertexPixels:
        vertexPixels.append(pixel)
        vertexCharge += pixVal
  return clusterCharge, vertexPixels, vertexCharge


def getMCPartE(ioll, tid):
  mctracks = ioll.get_data(larlite.data.kMCTrack, "mcreco")
  mcshowers = ioll.get_data(larlite.data.kMCShower, "mcreco")
  for mcparticles in [mctracks, mcshowers]:
    for mcpart in mcparticles:
      if mcpart.TrackID() == tid:
        return mcpart.Start().E()
  return -1.


def getMCProngParticle(sparseimg_vv, mcpg, mcpm, adc_v, ioll):

  particleDict = {}
  trackDict = {}
  totalPixI = 0.

  for p in range(3):
    for pix in sparseimg_vv[p]:
      totalPixI += pix.val
      pixContents = mcpm.getPixContent(p, pix.rawRow, pix.rawCol)
      for part in pixContents.particles:
        if abs(part.pdg) in particleDict:
          particleDict[abs(part.pdg)] += pixContents.pixI
        else:
          particleDict[abs(part.pdg)] = pixContents.pixI
        if part.tid in trackDict:
          trackDict[part.tid][2] += pixContents.pixI
        else:
          trackDict[part.tid] = [part.pdg, part.nodeidx, pixContents.pixI]

  maxPartPDG = 0 
  maxPartNID = -1
  maxPartTID = -1
  maxPartI = 0.
  maxPartComp = 0.
  maxPartE = -1.
  pdglist = []
  puritylist = []

  for part in particleDict:
    pdglist.append(part)
    puritylist.append(particleDict[part]/totalPixI)

  for track in trackDict:
    if trackDict[track][2] > maxPartI:
      maxPartI = trackDict[track][2]
      maxPartPDG = trackDict[track][0]
      maxPartNID = trackDict[track][1]
      maxPartTID = track

  totNodePixI = 0.
  if maxPartI > 0.:
    maxPartE = getMCPartE(ioll, maxPartTID)
    maxPartNode = mcpg.node_v[maxPartNID]
    if maxPartNode.tid != maxPartTID:
      sys.exit("ERROR: mismatch between node track id from mcpm and mcpg in getMCProngParticle")
    for p in range(3):
      pixels = maxPartNode.pix_vv[p]
      for iP in range(pixels.size()//2):
        row = (pixels[2*iP] - 2400)//6
        col = pixels[2*iP+1]
        totNodePixI += adc_v[p].pixel(row, col)
    if totNodePixI > 0.:
      maxPartComp = maxPartI/totNodePixI

  if maxPartComp > 1.:
    print("ERROR: prong completeness calculated to be >1")

  #return maxPartPDG, maxPartTID, totNodePixI, maxPartI/totalPixI, maxPartComp, pdglist, puritylist
  return maxPartPDG, maxPartTID, maxPartE, maxPartI/totalPixI, maxPartComp, pdglist, puritylist


def makeImage(prong_vv):
  plane0pix_row = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_col = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_val = np.zeros(prong_vv[0].size(), dtype=float)
  plane1pix_row = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_col = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_val = np.zeros(prong_vv[1].size(), dtype=float)
  plane2pix_row = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_col = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_val = np.zeros(prong_vv[2].size(), dtype=float)
  raw_plane0pix_row = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_col = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_val = np.zeros(prong_vv[3].size(), dtype=float)
  raw_plane1pix_row = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_col = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_val = np.zeros(prong_vv[4].size(), dtype=float)
  raw_plane2pix_row = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_col = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_val = np.zeros(prong_vv[5].size(), dtype=float)
  for i, pix in enumerate(prong_vv[0]):
    plane0pix_row[i] = pix.row
    plane0pix_col[i] = pix.col
    plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[1]):
    plane1pix_row[i] = pix.row
    plane1pix_col[i] = pix.col
    plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[2]):
    plane2pix_row[i] = pix.row
    plane2pix_col[i] = pix.col
    plane2pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[3]):
    raw_plane0pix_row[i] = pix.row
    raw_plane0pix_col[i] = pix.col
    raw_plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[4]):
    raw_plane1pix_row[i] = pix.row
    raw_plane1pix_col[i] = pix.col
    raw_plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[5]):
    raw_plane2pix_row[i] = pix.row
    raw_plane2pix_col[i] = pix.col
    raw_plane2pix_val[i] = pix.val
  image = np.zeros((6,512,512))
  image[0, plane0pix_row, plane0pix_col] = plane0pix_val
  image[2, plane1pix_row, plane1pix_col] = plane1pix_val
  image[4, plane2pix_row, plane2pix_col] = plane2pix_val
  image[1, raw_plane0pix_row, raw_plane0pix_col] = raw_plane0pix_val
  image[3, raw_plane1pix_row, raw_plane1pix_col] = raw_plane1pix_val
  image[5, raw_plane2pix_row, raw_plane2pix_col] = raw_plane2pix_val
  image = torch.from_numpy(image).float()
  norm = transforms.Normalize(mean, std)
  image = norm(image).reshape(1,6,512,512)
  return torch.clamp(image, max=4.0)


def getPID(cnnClass):
    if cnnClass == 0:
        return 11
    if cnnClass == 1:
        return 22
    if cnnClass == 2:
        return 13
    if cnnClass == 3:
        return 211
    if cnnClass == 4:
        return 2212
    return 0


#other needed classes
sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()
prongvars = larflow.reco.NuSelProngVars()
wcoverlapvars = larflow.reco.NuSelWCTaggerOverlap()
flowTriples = larflow.prep.FlowTriples()
piKEestimator = pionRange2T()
clusterFuncs = larflow.reco.ClusterFunctions()

model = ResNet34(2, ResBlock, outputs=5)
if "cuda" in args.device and args.multiGPU:
  model = nn.DataParallel(model)
if args.device == "cpu":
  checkpoint = torch.load(args.model_path, map_location=torch.device('cpu'))
else:
  checkpoint = torch.load(args.model_path)
try:
  model.load_state_dict(checkpoint['model_state_dict'])
except:
  model.module.load_state_dict(checkpoint['model_state_dict'])
model.to(args.device)
model.eval()

outRootFile = rt.TFile(args.outfile, "RECREATE")

if args.isMC:
  potTree = rt.TTree("potTree","potTree")
  totPOT = array('f', [0.])
  totGoodPOT = array('f', [0.])
  potTree.Branch("totPOT", totPOT, 'totPOT/F')
  potTree.Branch("totGoodPOT", totGoodPOT, 'totGoodPOT/F')

eventTree = rt.TTree("EventTree","EventTree")
maxNKpts = 500
maxNTrks = 100
maxNShwrs = 100
maxNParts = 1000
fileid = array('i', [0])
run = array('i', [0])
subrun = array('i', [0])
event = array('i', [0])
if args.isMC:
  xsecWeight = array('f', [0.])
  trueNuE = array('f', [0.])
  trueNuPDG = array('i', [0])
  trueNuCCNC = array('i', [0])
  trueNuMode = array('i', [0])
  trueNuIntrxnType = array('i', [0])
  trueVtxX = array('f', [0.])
  trueVtxY = array('f', [0.])
  trueVtxZ = array('f', [0.])
  trueLepE = array('f', [0.])
  trueLepPDG = array('i', [0])
  nTruePrimParts = array('i', [0])
  truePrimPartPDG = array('i', maxNParts*[0])
  truePrimPartX = array('f', maxNParts*[0.])
  truePrimPartY = array('f', maxNParts*[0.])
  truePrimPartZ = array('f', maxNParts*[0.])
  truePrimPartPx = array('f', maxNParts*[0.])
  truePrimPartPy = array('f', maxNParts*[0.])
  truePrimPartPz = array('f', maxNParts*[0.])
  truePrimPartE = array('f', maxNParts*[0.])
  truePrimPartContained = array('i', maxNParts*[0])
  nTrueSimParts = array('i', [0])
  trueSimPartPDG = array('i', maxNParts*[0])
  trueSimPartTID = array('i', maxNParts*[0])
  trueSimPartMID = array('i', maxNParts*[0])
  trueSimPartProcess = array('i', maxNParts*[0])
  trueSimPartX = array('f', maxNParts*[0.])
  trueSimPartY = array('f', maxNParts*[0.])
  trueSimPartZ = array('f', maxNParts*[0.])
  trueSimPartEDepX = array('f', maxNParts*[0.])
  trueSimPartEDepY = array('f', maxNParts*[0.])
  trueSimPartEDepZ = array('f', maxNParts*[0.])
  trueSimPartPx = array('f', maxNParts*[0.])
  trueSimPartPy = array('f', maxNParts*[0.])
  trueSimPartPz = array('f', maxNParts*[0.])
  trueSimPartE = array('f', maxNParts*[0.])
  trueSimPartEndPx = array('f', maxNParts*[0.])
  trueSimPartEndPy = array('f', maxNParts*[0.])
  trueSimPartEndPz = array('f', maxNParts*[0.])
  trueSimPartEndE = array('f', maxNParts*[0.])
  trueSimPartEndX = array('f', maxNParts*[0.])
  trueSimPartEndY = array('f', maxNParts*[0.])
  trueSimPartEndZ = array('f', maxNParts*[0.])
  trueSimPartContained = array('i', maxNParts*[0])
recoNuE = array('f', [0.])
foundVertex = array('i', [0])
vtxX = array('f', [0.])
vtxY = array('f', [0.])
vtxZ = array('f', [0.])
vtxIsFiducial = array('i', [0])
vtxContainment = array('i', [0])
if args.isMC:
  vtxDistToTrue = array('f', [0.])
vtxScore = array('f', [0.])
vtxFracHitsOnCosmic = array('f', [0.])
eventPCAxis0 = array('f', 3*[0.])
eventPCAxis1 = array('f', 3*[0.])
eventPCAxis2 = array('f', 3*[0.])
eventPCAxis0TSlope = array('i', [0])
eventPCEigenVals = array('f', 3*[0.])
eventPCProjMaxGap = array('f', 5*[0.])
eventPCProjMaxDist = array('f', 5*[0.])
if not args.noKeypoints:
  nKeypoints = array('i', [0])
  kpClusterType = array('i', maxNKpts*[0])
  kpFilterType = array('i', maxNKpts*[0])
  kpMaxScore = array('f', maxNKpts*[0.])
  kpMaxPosX = array('f', maxNKpts*[0.])
  kpMaxPosY = array('f', maxNKpts*[0.])
  kpMaxPosZ = array('f', maxNKpts*[0.])
nTracks = array('i', [0])
trackIsSecondary = array('i', maxNTrks*[0])
trackNHits = array('i', maxNTrks*[0])
trackHitFrac = array('f', maxNTrks*[0.])
trackCharge = array('f', maxNTrks*[0.])
trackChargeFrac = array('f', maxNTrks*[0.])
trackCosTheta = array('f', maxNTrks*[0.])
trackCosThetaY = array('f', maxNTrks*[0.])
trackDistToVtx = array('f', maxNTrks*[0.])
trackStartPosX = array('f', maxNTrks*[0.])
trackStartPosY = array('f', maxNTrks*[0.])
trackStartPosZ = array('f', maxNTrks*[0.])
trackStartDirX = array('f', maxNTrks*[0.])
trackStartDirY = array('f', maxNTrks*[0.])
trackStartDirZ = array('f', maxNTrks*[0.])
trackEndPosX = array('f', maxNTrks*[0.])
trackEndPosY = array('f', maxNTrks*[0.])
trackEndPosZ = array('f', maxNTrks*[0.])
trackClassified = array('i', maxNTrks*[0])
trackPID = array('i', maxNTrks*[0])
trackElScore = array('f', maxNTrks*[0.])
trackPhScore = array('f', maxNTrks*[0.])
trackMuScore = array('f', maxNTrks*[0.])
trackPiScore = array('f', maxNTrks*[0.])
trackPrScore = array('f', maxNTrks*[0.])
trackComp = array('f', maxNTrks*[0.])
trackPurity = array('f', maxNTrks*[0.])
trackProcess = array('i', maxNTrks*[0])
trackPrimaryScore = array('f', maxNTrks*[0.])
trackFromNeutralScore = array('f', maxNTrks*[0.])
trackFromChargedScore = array('f', maxNTrks*[0.])
trackPrimScore = array('f', maxNTrks*[0.])
trackRecoE = array('f', maxNTrks*[0.])
if args.isMC:
  trackTruePID = array('i', maxNTrks*[0])
  trackTrueTID = array('i', maxNTrks*[0])
  trackTrueE = array('f', maxNTrks*[0])
  trackTruePurity = array('f', maxNTrks*[0.])
  trackTrueComp = array('f', maxNTrks*[0.])
  trackTrueElPurity = array('f', maxNTrks*[0.])
  trackTruePhPurity = array('f', maxNTrks*[0.])
  trackTrueMuPurity = array('f', maxNTrks*[0.])
  trackTruePiPurity = array('f', maxNTrks*[0.])
  trackTruePrPurity = array('f', maxNTrks*[0.])
nShowers = array('i', [0])
showerIsSecondary = array('i', maxNShwrs*[0])
showerNHits = array('i', maxNShwrs*[0])
showerHitFrac = array('f', maxNShwrs*[0.])
showerCharge = array('f', maxNShwrs*[0.])
showerChargeFrac = array('f', maxNShwrs*[0.])
showerCosTheta = array('f', maxNShwrs*[0.])
showerCosThetaY = array('f', maxNShwrs*[0.])
showerDistToVtx = array('f', maxNShwrs*[0.])
showerStartPosX = array('f', maxNShwrs*[0.])
showerStartPosY = array('f', maxNShwrs*[0.])
showerStartPosZ = array('f', maxNShwrs*[0.])
showerStartDirX = array('f', maxNShwrs*[0.])
showerStartDirY = array('f', maxNShwrs*[0.])
showerStartDirZ = array('f', maxNShwrs*[0.])
showerClassified = array('i', maxNShwrs*[0])
showerPID = array('i', maxNShwrs*[0])
showerElScore = array('f', maxNShwrs*[0.])
showerPhScore = array('f', maxNShwrs*[0.])
showerMuScore = array('f', maxNShwrs*[0.])
showerPiScore = array('f', maxNShwrs*[0.])
showerPrScore = array('f', maxNShwrs*[0.])
showerComp = array('f', maxNShwrs*[0])
showerPurity = array('f', maxNShwrs*[0])
showerProcess = array('i', maxNShwrs*[0])
showerPrimaryScore = array('f', maxNShwrs*[0.])
showerFromNeutralScore = array('f', maxNShwrs*[0.])
showerFromChargedScore = array('f', maxNShwrs*[0.])
showerRecoE = array('f', maxNShwrs*[0])
if args.isMC:
  showerTruePID = array('i', maxNShwrs*[0])
  showerTrueTID = array('i', maxNShwrs*[0])
  showerTrueE = array('f', maxNShwrs*[0])
  showerTruePurity = array('f', maxNShwrs*[0.])
  showerTrueComp = array('f', maxNShwrs*[0.])
  showerTrueElPurity = array('f', maxNShwrs*[0.])
  showerTruePhPurity = array('f', maxNShwrs*[0.])
  showerTrueMuPurity = array('f', maxNShwrs*[0.])
  showerTruePiPurity = array('f', maxNShwrs*[0.])
  showerTruePrPurity = array('f', maxNShwrs*[0.])
eventTree.Branch("fileid", fileid, 'fileid/I')
eventTree.Branch("run", run, 'run/I')
eventTree.Branch("subrun", subrun, 'subrun/I')
eventTree.Branch("event", event, 'event/I')
if args.isMC:
  eventTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
  eventTree.Branch("trueNuE", trueNuE, 'trueNuE/F')
  eventTree.Branch("trueNuPDG", trueNuPDG, 'trueNuPDG/I')
  eventTree.Branch("trueNuCCNC", trueNuCCNC, 'trueNuCCNC/I')
  eventTree.Branch("trueNuMode", trueNuMode, 'trueNuMode/I')
  eventTree.Branch("trueNuIntrxnType", trueNuIntrxnType, 'trueNuIntrxnType/I')
  eventTree.Branch("trueVtxX", trueVtxX, 'trueVtxX/F')
  eventTree.Branch("trueVtxY", trueVtxY, 'trueVtxY/F')
  eventTree.Branch("trueVtxZ", trueVtxZ, 'trueVtxZ/F')
  eventTree.Branch("trueLepE", trueLepE, 'trueLepE/F')
  eventTree.Branch("trueLepPDG", trueLepPDG, 'trueLepPDG/I')
  eventTree.Branch("nTruePrimParts", nTruePrimParts, 'nTruePrimParts/I')
  eventTree.Branch("truePrimPartPDG", truePrimPartPDG, 'truePrimPartPDG[nTruePrimParts]/I')
  eventTree.Branch("truePrimPartX", truePrimPartX, 'truePrimPartX[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartY", truePrimPartY, 'truePrimPartY[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartZ", truePrimPartZ, 'truePrimPartZ[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartPx", truePrimPartPx, 'truePrimPartPx[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartPy", truePrimPartPy, 'truePrimPartPy[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartPz", truePrimPartPz, 'truePrimPartPz[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartE", truePrimPartE, 'truePrimPartE[nTruePrimParts]/F')
  eventTree.Branch("truePrimPartContained", truePrimPartContained, 'truePrimPartContained[nTruePrimParts]/I')
  eventTree.Branch("nTrueSimParts", nTrueSimParts, 'nTrueSimParts/I')
  eventTree.Branch("trueSimPartPDG", trueSimPartPDG, 'trueSimPartPDG[nTrueSimParts]/I')
  eventTree.Branch("trueSimPartTID", trueSimPartTID, 'trueSimPartTID[nTrueSimParts]/I')
  eventTree.Branch("trueSimPartMID", trueSimPartMID, 'trueSimPartMID[nTrueSimParts]/I')
  eventTree.Branch("trueSimPartProcess", trueSimPartProcess, 'trueSimPartProcess[nTrueSimParts]/I')
  eventTree.Branch("trueSimPartX", trueSimPartX, 'trueSimPartX[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartY", trueSimPartY, 'trueSimPartY[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartZ", trueSimPartZ, 'trueSimPartZ[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEDepX", trueSimPartEDepX, 'trueSimPartEDepX[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEDepY", trueSimPartEDepY, 'trueSimPartEDepY[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEDepZ", trueSimPartEDepZ, 'trueSimPartEDepZ[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartPx", trueSimPartPx, 'trueSimPartPx[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartPy", trueSimPartPy, 'trueSimPartPy[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartPz", trueSimPartPz, 'trueSimPartPz[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartE", trueSimPartE, 'trueSimPartE[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndPx", trueSimPartEndPx, 'trueSimPartEndPx[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndPy", trueSimPartEndPy, 'trueSimPartEndPy[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndPz", trueSimPartEndPz, 'trueSimPartEndPz[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndE", trueSimPartEndE, 'trueSimPartEndE[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndX", trueSimPartEndX, 'trueSimPartEndX[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndY", trueSimPartEndY, 'trueSimPartEndY[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartEndZ", trueSimPartEndZ, 'trueSimPartEndZ[nTrueSimParts]/F')
  eventTree.Branch("trueSimPartContained", trueSimPartContained, 'trueSimPartContained[nTrueSimParts]/I')
eventTree.Branch("recoNuE", recoNuE, 'recoNuE/F')
eventTree.Branch("foundVertex", foundVertex, 'foundVertex/I')
eventTree.Branch("vtxX", vtxX, 'vtxX/F')
eventTree.Branch("vtxY", vtxY, 'vtxY/F')
eventTree.Branch("vtxZ", vtxZ, 'vtxZ/F')
eventTree.Branch("vtxIsFiducial", vtxIsFiducial, 'vtxIsFiducial/I')
eventTree.Branch("vtxContainment", vtxContainment, 'vtxContainment/I')
if args.isMC:
  eventTree.Branch("vtxDistToTrue", vtxDistToTrue, 'vtxDistToTrue/F')
eventTree.Branch("vtxScore", vtxScore, 'vtxScore/F')
eventTree.Branch("vtxFracHitsOnCosmic", vtxFracHitsOnCosmic, 'vtxFracHitsOnCosmic/F')
eventTree.Branch("eventPCAxis0", eventPCAxis0, 'eventPCAxis0[3]/F')
eventTree.Branch("eventPCAxis1", eventPCAxis1, 'eventPCAxis1[3]/F')
eventTree.Branch("eventPCAxis2", eventPCAxis2, 'eventPCAxis2[3]/F')
eventTree.Branch("eventPCAxis0TSlope", eventPCAxis0TSlope, 'eventPCAxis0TSlope/I')
eventTree.Branch("eventPCEigenVals", eventPCEigenVals, 'eventPCEigenVals[3]/F')
eventTree.Branch("eventPCProjMaxGap", eventPCProjMaxGap, 'eventPCProjMaxGap[5]/F')
eventTree.Branch("eventPCProjMaxDist", eventPCProjMaxDist, 'eventPCProjMaxDist[5]/F')
if not args.noKeypoints:
  eventTree.Branch("nKeypoints", nKeypoints, 'nKeypoints/I')
  eventTree.Branch("kpClusterType", kpClusterType, 'kpClusterType[nKeypoints]/I')
  eventTree.Branch("kpFilterType", kpFilterType, 'kpFilterType[nKeypoints]/I')
  eventTree.Branch("kpMaxScore", kpMaxScore, 'kpMaxScore[nKeypoints]/F')
  eventTree.Branch("kpMaxPosX", kpMaxPosX, 'kpMaxPosX[nKeypoints]/F')
  eventTree.Branch("kpMaxPosY", kpMaxPosY, 'kpMaxPosY[nKeypoints]/F')
  eventTree.Branch("kpMaxPosZ", kpMaxPosZ, 'kpMaxPosZ[nKeypoints]/F')
eventTree.Branch("nTracks", nTracks, 'nTracks/I')
eventTree.Branch("trackIsSecondary", trackIsSecondary, 'trackIsSecondary[nTracks]/I')
eventTree.Branch("trackNHits", trackNHits, 'trackNHits[nTracks]/I')
eventTree.Branch("trackHitFrac", trackHitFrac, 'trackHitFrac[nTracks]/F')
eventTree.Branch("trackCharge", trackCharge, 'trackCharge[nTracks]/F')
eventTree.Branch("trackChargeFrac", trackChargeFrac, 'trackChargeFrac[nTracks]/F')
eventTree.Branch("trackCosTheta", trackCosTheta, 'trackCosTheta[nTracks]/F')
eventTree.Branch("trackCosThetaY", trackCosThetaY, 'trackCosThetaY[nTracks]/F')
eventTree.Branch("trackDistToVtx", trackDistToVtx, 'trackDistToVtx[nTracks]/F')
eventTree.Branch("trackStartPosX", trackStartPosX, 'trackStartPosX[nTracks]/F')
eventTree.Branch("trackStartPosY", trackStartPosY, 'trackStartPosY[nTracks]/F')
eventTree.Branch("trackStartPosZ", trackStartPosZ, 'trackStartPosZ[nTracks]/F')
eventTree.Branch("trackStartDirX", trackStartDirX, 'trackStartDirX[nTracks]/F')
eventTree.Branch("trackStartDirY", trackStartDirY, 'trackStartDirY[nTracks]/F')
eventTree.Branch("trackStartDirZ", trackStartDirZ, 'trackStartDirZ[nTracks]/F')
eventTree.Branch("trackEndPosX", trackEndPosX, 'trackEndPosX[nTracks]/F')
eventTree.Branch("trackEndPosY", trackEndPosY, 'trackEndPosY[nTracks]/F')
eventTree.Branch("trackEndPosZ", trackEndPosZ, 'trackEndPosZ[nTracks]/F')
eventTree.Branch("trackClassified", trackClassified, 'trackClassified[nTracks]/I')
eventTree.Branch("trackPID", trackPID, 'trackPID[nTracks]/I')
eventTree.Branch("trackElScore", trackElScore, 'trackElScore[nTracks]/F')
eventTree.Branch("trackPhScore", trackPhScore, 'trackPhScore[nTracks]/F')
eventTree.Branch("trackMuScore", trackMuScore, 'trackMuScore[nTracks]/F')
eventTree.Branch("trackPiScore", trackPiScore, 'trackPiScore[nTracks]/F')
eventTree.Branch("trackPrScore", trackPrScore, 'trackPrScore[nTracks]/F')
eventTree.Branch("trackComp", trackComp, 'trackComp[nTracks]/F')
eventTree.Branch("trackPurity", trackPurity, 'trackPurity[nTracks]/F')
eventTree.Branch("trackProcess", trackProcess, 'trackProcess[nTracks]/I')
eventTree.Branch("trackPrimaryScore", trackPrimaryScore, 'trackPrimaryScore[nTracks]/F')
eventTree.Branch("trackFromNeutralScore", trackFromNeutralScore, 'trackFromNeutralScore[nTracks]/F')
eventTree.Branch("trackFromChargedScore", trackFromChargedScore, 'trackFromChargedScore[nTracks]/F')
eventTree.Branch("trackRecoE", trackRecoE, 'trackRecoE[nTracks]/F')
if args.isMC:
  eventTree.Branch("trackTruePID", trackTruePID, 'trackTruePID[nTracks]/I')
  eventTree.Branch("trackTrueTID", trackTrueTID, 'trackTrueTID[nTracks]/I')
  eventTree.Branch("trackTrueE", trackTrueE, 'trackTrueE[nTracks]/F')
  eventTree.Branch("trackTruePurity", trackTruePurity, 'trackTruePurity[nTracks]/F')
  eventTree.Branch("trackTrueComp", trackTrueComp, 'trackTrueComp[nTracks]/F')
  eventTree.Branch("trackTrueElPurity", trackTrueElPurity, 'trackTrueElPurity[nTracks]/F')
  eventTree.Branch("trackTruePhPurity", trackTruePhPurity, 'trackTruePhPurity[nTracks]/F')
  eventTree.Branch("trackTrueMuPurity", trackTrueMuPurity, 'trackTrueMuPurity[nTracks]/F')
  eventTree.Branch("trackTruePiPurity", trackTruePiPurity, 'trackTruePiPurity[nTracks]/F')
  eventTree.Branch("trackTruePrPurity", trackTruePrPurity, 'trackTruePrPurity[nTracks]/F')
eventTree.Branch("nShowers", nShowers, 'nShowers/I')
eventTree.Branch("showerIsSecondary", showerIsSecondary, 'showerIsSecondary[nShowers]/I')
eventTree.Branch("showerNHits", showerNHits, 'showerNHits[nShowers]/I')
eventTree.Branch("showerHitFrac", showerHitFrac, 'showerHitFrac[nShowers]/F')
eventTree.Branch("showerCharge", showerCharge, 'showerCharge[nShowers]/F')
eventTree.Branch("showerChargeFrac", showerChargeFrac, 'showerChargeFrac[nShowers]/F')
eventTree.Branch("showerCosTheta", showerCosTheta, 'showerCosTheta[nShowers]/F')
eventTree.Branch("showerCosThetaY", showerCosThetaY, 'showerCosThetaY[nShowers]/F')
eventTree.Branch("showerDistToVtx", showerDistToVtx, 'showerDistToVtx[nShowers]/F')
eventTree.Branch("showerStartPosX", showerStartPosX, 'showerStartPosX[nShowers]/F')
eventTree.Branch("showerStartPosY", showerStartPosY, 'showerStartPosY[nShowers]/F')
eventTree.Branch("showerStartPosZ", showerStartPosZ, 'showerStartPosZ[nShowers]/F')
eventTree.Branch("showerStartDirX", showerStartDirX, 'showerStartDirX[nShowers]/F')
eventTree.Branch("showerStartDirY", showerStartDirY, 'showerStartDirY[nShowers]/F')
eventTree.Branch("showerStartDirZ", showerStartDirZ, 'showerStartDirZ[nShowers]/F')
eventTree.Branch("showerClassified", showerClassified, 'showerClassified[nShowers]/I')
eventTree.Branch("showerPID", showerPID, 'showerPID[nShowers]/I')
eventTree.Branch("showerElScore", showerElScore, 'showerElScore[nShowers]/F')
eventTree.Branch("showerPhScore", showerPhScore, 'showerPhScore[nShowers]/F')
eventTree.Branch("showerMuScore", showerMuScore, 'showerMuScore[nShowers]/F')
eventTree.Branch("showerPiScore", showerPiScore, 'showerPiScore[nShowers]/F')
eventTree.Branch("showerPrScore", showerPrScore, 'showerPrScore[nShowers]/F')
eventTree.Branch("showerComp", showerComp, 'showerComp[nShowers]/F')
eventTree.Branch("showerPurity", showerPurity, 'showerPurity[nShowers]/F')
eventTree.Branch("showerProcess", showerProcess, 'showerProcess[nShowers]/I')
eventTree.Branch("showerPrimaryScore", showerPrimaryScore, 'showerPrimaryScore[nShowers]/F')
eventTree.Branch("showerFromNeutralScore", showerFromNeutralScore, 'showerFromNeutralScore[nShowers]/F')
eventTree.Branch("showerFromChargedScore", showerFromChargedScore, 'showerFromChargedScore[nShowers]/F')
eventTree.Branch("showerRecoE", showerRecoE, 'showerRecoE[nShowers]/F')
if args.isMC:
  eventTree.Branch("showerTruePID", showerTruePID, 'showerTruePID[nShowers]/I')
  eventTree.Branch("showerTrueTID", showerTrueTID, 'showerTrueTID[nShowers]/I')
  eventTree.Branch("showerTrueE", showerTrueE, 'showerTrueE[nShowers]/F')
  eventTree.Branch("showerTruePurity", showerTruePurity, 'showerTruePurity[nShowers]/F')
  eventTree.Branch("showerTrueComp", showerTrueComp, 'showerTrueComp[nShowers]/F')
  eventTree.Branch("showerTrueElPurity", showerTrueElPurity, 'showerTrueElPurity[nShowers]/F')
  eventTree.Branch("showerTruePhPurity", showerTruePhPurity, 'showerTruePhPurity[nShowers]/F')
  eventTree.Branch("showerTrueMuPurity", showerTrueMuPurity, 'showerTrueMuPurity[nShowers]/F')
  eventTree.Branch("showerTruePiPurity", showerTruePiPurity, 'showerTruePiPurity[nShowers]/F')
  eventTree.Branch("showerTruePrPurity", showerTruePrPurity, 'showerTruePrPurity[nShowers]/F')


if args.isMC:
  totPOT_ = 0.
  totGoodPOT_ = 0.
  if not args.ignoreWeights:
    weights = Weights(args.weightfile)


#-------- begin file loop -----------------------------------------------------#
for filepair in files:

  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(filepair[1])
  ioll.open()

  iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv.add_in_file(filepair[1])
  iolcv.reverse_all_products()
  iolcv.initialize()

  kpsfile = rt.TFile(filepair[0])
  kpst = kpsfile.Get("KPSRecoManagerTree")

  try:
    nKPSTEntries = kpst.GetEntries()
  except:
    print("WARNING: %s is empty. skipping..."%(filepair[0]))
    ioll.close()
    iolcv.finalize()
    kpsfile.Close()
    continue

  if args.isMC:
    potInFile, goodPotInFile = SumPOT(filepair[1])
    if potInFile < 0. and goodPotInFile < 0.:
      print("WARNING: merged dlreco/dlana file does not have POT info. Skipping this input")
      continue
    totPOT_ = totPOT_ + potInFile
    totGoodPOT_ = totGoodPOT_ + goodPotInFile

  #++++++ begin entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
  for ientry in range(ioll.get_entries()):

    #print("reached entry:", ientry)

    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    kpst.GetEntry(ientry)
  
    if kpst.run != ioll.run_id() or kpst.subrun != ioll.subrun_id() or kpst.event != ioll.event_id():
      print("WARNING: EVENTS DON'T MATCH!!!")
      print("truth run/subrun/event: %i/%i/%i"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
      print("reco run/subrun/event: %i/%i/%i"%(kpst.run,kpst.subrun,kpst.event))
      continue


    if args.isMC:

      mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
      nuInt = mctruth.at(0).GetNeutrino()
      lep = nuInt.Lepton()
      mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
      trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

      if not isFiducialWCSCE(trueVtxPos):
        continue

      if args.ignoreWeights:
        xsecWeight[0] = 1.
      else:
        try:
          xsecWeight[0] = weights.get(kpst.run, kpst.subrun, kpst.event)
          if isinf(xsecWeight[0]):
            continue
        except:
          if args.skipNoWeightEvts:
            print("WARNING: Couldn't find xsec weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))
            continue
          else:
            sys.exit("ERROR: Couldn't find xsec weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))

      if nuInt.CCNC() == 0:
        trueLepPDG[0] = lep.PdgCode()
        trueLepE[0] = lep.Momentum().E()
      else:
        trueLepPDG[0] = 0
        trueLepE[0] = -9.

      trueNuPDG[0] = nuInt.Nu().PdgCode()
      trueNuCCNC[0] = nuInt.CCNC()
      trueNuMode[0] = nuInt.Mode()
      trueNuIntrxnType[0] = nuInt.InteractionType()
      trueNuE[0] = nuInt.Nu().Momentum().E()
      trueVtxX[0] = trueVtxPos.X()
      trueVtxY[0] = trueVtxPos.Y()
      trueVtxZ[0] = trueVtxPos.Z()

      nTruePrimParts[0] = 0
      iPP = 0
      for mcpart in mctruth.at(0).GetParticles():
        if mcpart.StatusCode() == 1:
          nTruePrimParts[0] += 1
          truePrimPartPDG[iPP] = mcpart.PdgCode()
          sceCorrectedStartPos = getSCECorrectedPos(mcpart.Position(0), sce)
          sceCorrectedEndPos = getSCECorrectedPos(mcpart.Position(mcpart.Trajectory().size()-1), sce)
          truePrimPartX[iPP] = sceCorrectedStartPos.X()
          truePrimPartY[iPP] = sceCorrectedStartPos.Y()
          truePrimPartZ[iPP] = sceCorrectedStartPos.Z()
          truePrimPartPx[iPP] = mcpart.Momentum(0).Px()
          truePrimPartPy[iPP] = mcpart.Momentum(0).Py()
          truePrimPartPz[iPP] = mcpart.Momentum(0).Pz()
          truePrimPartE[iPP] = mcpart.Momentum(0).E()
          truePrimPartContained[iPP] = isFiducialWCSCE(sceCorrectedEndPos)
          iPP += 1

      mctracks = ioll.get_data(larlite.data.kMCTrack, "mcreco")
      mcshowers = ioll.get_data(larlite.data.kMCShower, "mcreco")

      nTrueSimParts[0] = 0
      iDS = 0
      for mcparts in [mctracks, mcshowers]:
        for mcpart in mcparts:
          nTrueSimParts[0] += 1
          trueSimPartPDG[iDS] = mcpart.PdgCode()
          trueSimPartTID[iDS] = mcpart.TrackID()
          trueSimPartMID[iDS] = mcpart.MotherTrackID()
          if mcpart.Process() == 'primary':
            trueSimPartProcess[iDS] = 0
          elif mcpart.Process() == 'Decay':
            trueSimPartProcess[iDS] = 1
          else:
            trueSimPartProcess[iDS] = 2
          sceCorrectedStartPos = getSCECorrectedPos(mcpart.Start(), sce)
          sceCorrectedEndPos = getSCECorrectedPos(mcpart.End(), sce)
          trueSimPartX[iDS] = sceCorrectedStartPos.X()
          trueSimPartY[iDS] = sceCorrectedStartPos.Y()
          trueSimPartZ[iDS] = sceCorrectedStartPos.Z()
          if mcpart.PdgCode() == 22:
            trueSimPartEDepX[iDS] = mcpart.DetProfile().X()
            trueSimPartEDepY[iDS] = mcpart.DetProfile().Y()
            trueSimPartEDepZ[iDS] = mcpart.DetProfile().Z()
          else:
            trueSimPartEDepX[iDS] = sceCorrectedStartPos.X()
            trueSimPartEDepY[iDS] = sceCorrectedStartPos.Y()
            trueSimPartEDepZ[iDS] = sceCorrectedStartPos.Z()
          trueSimPartPx[iDS] = mcpart.Start().Px()
          trueSimPartPy[iDS] = mcpart.Start().Py()
          trueSimPartPz[iDS] = mcpart.Start().Pz()
          trueSimPartE[iDS] = mcpart.Start().E()
          trueSimPartEndPx[iDS] = mcpart.End().Px()
          trueSimPartEndPy[iDS] = mcpart.End().Py()
          trueSimPartEndPz[iDS] = mcpart.End().Pz()
          trueSimPartEndE[iDS] = mcpart.End().E()
          trueSimPartEndX[iDS] = sceCorrectedEndPos.X()
          trueSimPartEndY[iDS] = sceCorrectedEndPos.Y()
          trueSimPartEndZ[iDS] = sceCorrectedEndPos.Z()
          trueSimPartContained[iDS] = isFiducialWCSCE(sceCorrectedEndPos)
          iDS += 1

    #else: #from "if args.isMC"
    #  xsecWeight[0] = -1.
    #  trueLepPDG[0] = 0
    #  trueLepE[0] = -9.
    #  trueNuPDG[0] = 0
    #  trueNuCCNC[0] = -1
    #  trueNuE[0] = -9.
    #  trueVtxX[0] = -999.
    #  trueVtxY[0] = -999.
    #  trueVtxZ[0] = -999.
    #  nTruePrimParts[0] = -1
    #  nTrueSimParts[0] = -1

    fileid[0] = -1
    for tag in filepair[0].split("_"):
      if 'fileid' in tag:
        fileid[0] = int(tag.replace("fileid",""))
        break
    run[0] = kpst.run
    subrun[0] = kpst.subrun
    event[0] = kpst.event

    if not args.noKeypoints:
      nKeypoints[0] = 0
      for kp in kpst.kpc_nu_v:
        kpClusterType[nKeypoints[0]] = kp._cluster_type
        kpFilterType[nKeypoints[0]] = 0
        kpMaxScore[nKeypoints[0]] = kp.max_score
        kpMaxPosX[nKeypoints[0]] = kp.max_pt_v[0]
        kpMaxPosY[nKeypoints[0]] = kp.max_pt_v[1]
        kpMaxPosZ[nKeypoints[0]] = kp.max_pt_v[2]
        nKeypoints[0] += 1
      for kp in kpst.kpc_track_v:
        kpClusterType[nKeypoints[0]] = kp._cluster_type
        kpFilterType[nKeypoints[0]] = 0
        kpMaxScore[nKeypoints[0]] = kp.max_score
        kpMaxPosX[nKeypoints[0]] = kp.max_pt_v[0]
        kpMaxPosY[nKeypoints[0]] = kp.max_pt_v[1]
        kpMaxPosZ[nKeypoints[0]] = kp.max_pt_v[2]
        nKeypoints[0] += 1
      for kp in kpst.kpc_shower_v:
        kpClusterType[nKeypoints[0]] = kp._cluster_type
        kpFilterType[nKeypoints[0]] = 0
        kpMaxScore[nKeypoints[0]] = kp.max_score
        kpMaxPosX[nKeypoints[0]] = kp.max_pt_v[0]
        kpMaxPosY[nKeypoints[0]] = kp.max_pt_v[1]
        kpMaxPosZ[nKeypoints[0]] = kp.max_pt_v[2]
        nKeypoints[0] += 1
      for kp in kpst.kpc_cosmic_v:
        kpClusterType[nKeypoints[0]] = kp._cluster_type
        kpFilterType[nKeypoints[0]] = 1
        kpMaxScore[nKeypoints[0]] = kp.max_score
        kpMaxPosX[nKeypoints[0]] = kp.max_pt_v[0]
        kpMaxPosY[nKeypoints[0]] = kp.max_pt_v[1]
        kpMaxPosZ[nKeypoints[0]] = kp.max_pt_v[2]
        nKeypoints[0] += 1

    foundVertex[0] = 0
    vtxScore[0] = -1.
    for vtx in kpst.nuvetoed_v:
      if vtx.keypoint_type != 0:
        continue
      foundVertex[0] = 1
      if vtx.netNuScore > vtxScore[0]:
        vtxScore[0] = vtx.netNuScore
        vertex = vtx

    if foundVertex[0] == 0:
      recoNuE[0] = -9.
      vtxX[0] = -999.
      vtxY[0] = -999.
      vtxZ[0] = -999.
      vtxIsFiducial[0] = -1
      vtxContainment[0] = -1
      if args.isMC:
        vtxDistToTrue[0] = -99.
      vtxFracHitsOnCosmic[0] = -1.
      nTracks[0] = 0
      nShowers[0] = 0
      for iPCA in range(3):
        eventPCAxis0[iPCA] = 0.
        eventPCAxis1[iPCA] = 0.
        eventPCAxis2[iPCA] = 0.
        eventPCAxis0TSlope[0] = 0
        eventPCEigenVals[iPCA] = 0.
      for iPrj in range(5):
        eventPCProjMaxGap[iPrj] = -9.
        eventPCProjMaxDist[iPrj] = -9.
      eventTree.Fill()
      continue

    if args.isMC:
      vtxDistToTrue[0] = getVertexDistance(trueVtxPos, vertex)
      mcpg = ublarcvapp.mctools.MCPixelPGraph()
      mcpg.set_adc_treename("wire")
      mcpg.buildgraph(iolcv, ioll)
      mcpm = ublarcvapp.mctools.MCPixelPMap()
      mcpm.set_adc_treename("wire")
      mcpm.buildmap(iolcv, mcpg)
    #else:
    #  vtxDistToTrue[0] = -99.

    vtxX[0] = vertex.pos[0]
    vtxY[0] = vertex.pos[1]
    vtxZ[0] = vertex.pos[2]
    vtxTVec3 = rt.TVector3(vertex.pos[0], vertex.pos[1], vertex.pos[2])
    vtxIsFiducial[0] = int(isFiducialWCSCE(vtxTVec3))
    if vtxIsFiducial[0] == 0:
      vtxContainment[0] = 0
      prongsAreContained = False
    else:
      prongsAreContained = True

    nusel = larflow.reco.NuSelectionVariables()
    wcoverlapvars.analyze(vertex, nusel, iolcv)
    vtxFracHitsOnCosmic[0] = nusel.frac_allhits_on_cosmic

    nTracks[0] = vertex.track_v.size()
    nShowers[0] = vertex.shower_v.size()

    evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
    csmImage2D = iolcv.get_data(larcv.kProductImage2D, "thrumu")
    adc_v = evtImage2D.Image2DArray()
    thrumu_v = csmImage2D.Image2DArray()

    vertexPixels = []
    vertexCharge = 0.
    vertexNHits = 0

    eventLarflowCluster = larlite.larflowcluster()

    recoNuE[0] = 0.


    #++++++ begin track loop ++++++++++++++++++++++++++++++++++++++++++++++++++=
    for iTrk, trackCls in enumerate(vertex.track_hitcluster_v):

      for hit in trackCls:
        eventLarflowCluster.push_back(hit)
        if prongsAreContained:
          hitPt = rt.TVector3(hit[0], hit[1], hit[2])
          prongsAreContained = isFiducialWCSCE(hitPt)

      vertexNHits += trackCls.size()
      trackIsSecondary[iTrk] = vertex.track_isSecondary_v[iTrk]
      trackNHits[iTrk] = trackCls.size()
      trackCharge[iTrk], vertexPixels, vertexCharge = addClusterCharge(iolcv,trackCls,vertexPixels,vertexCharge,10.)
      nTrajPoints = vertex.track_v[iTrk].NumberTrajectoryPoints()
      trackLength = getDistance(vertex.track_v[iTrk].Vertex(),vertex.track_v[iTrk].End()) if (nTrajPoints > 1) else -9.
      goodTrack = nTrajPoints > 1 and trackLength > 1e-6
      trackCosTheta[iTrk] = getCosThetaBeamTrack(vertex.track_v[iTrk]) if goodTrack else -9.
      trackCosThetaY[iTrk] = getCosThetaGravTrack(vertex.track_v[iTrk]) if goodTrack else -9.
      trackDistToVtx[iTrk] = getVertexDistance(vertex.track_v[iTrk].Vertex(), vertex) if goodTrack else -9.
      trackStartPosX[iTrk] = vertex.track_v[iTrk].Vertex().X() if goodTrack else -9.
      trackStartPosY[iTrk] = vertex.track_v[iTrk].Vertex().Y() if goodTrack else -9.
      trackStartPosZ[iTrk] = vertex.track_v[iTrk].Vertex().Z() if goodTrack else -9.
      trackEndPosX[iTrk] = vertex.track_v[iTrk].End().X() if goodTrack else -9.
      trackEndPosY[iTrk] = vertex.track_v[iTrk].End().Y() if goodTrack else -9.
      trackEndPosZ[iTrk] = vertex.track_v[iTrk].End().Z() if goodTrack else -9.
      trackDirPt1 = vertex.track_v[iTrk].Vertex()
      trackDirPt2 = vertex.track_v[iTrk].Vertex()
      for iTrj in range(vertex.track_v[iTrk].NumberTrajectoryPoints()):
        trackDirPt2 = vertex.track_v[iTrk].LocationAtPoint(iTrj)
        if getDistance(trackDirPt1, trackDirPt2) > 5.:
          break
      trackDir = getDirection(trackDirPt1, trackDirPt2) if goodTrack else (0,0,0)
      trackStartDirX[iTrk], trackStartDirY[iTrk], trackStartDirZ[iTrk] = trackDir
      

      skip = True
      if goodTrack:
        skip = False
        cropPt = vertex.track_v[iTrk].End()
        prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,trackCls,cropPt,10.,512,512)
        for p in range(3):
          if prong_vv[p].size() < 10:
            skip = True
            break
      if skip:
        trackClassified[iTrk] = 0
        trackPID[iTrk] = 0
        trackElScore[iTrk] = -99.
        trackPhScore[iTrk] = -99.
        trackMuScore[iTrk] = -99.
        trackPiScore[iTrk] = -99.
        trackPrScore[iTrk] = -99.
        trackComp[iTrk] = -1.
        trackPurity[iTrk] = -1.
        trackProcess[iTrk] = -1
        trackPrimaryScore[iTrk] = -99.
        trackFromNeutralScore[iTrk] = -99.
        trackFromChargedScore[iTrk] = -99.
        trackRecoE[iTrk] = -1.
        continue

      with torch.no_grad():
        prongImage = makeImage(prong_vv).to(args.device)
        prongCNN_out = model(prongImage)
      trackClassified[iTrk] = 1
      trackPID[iTrk] = getPID(prongCNN_out[0].argmax(1).item())
      trackElScore[iTrk] = prongCNN_out[0][0][0].item()
      trackPhScore[iTrk] = prongCNN_out[0][0][1].item()
      trackMuScore[iTrk] = prongCNN_out[0][0][2].item()
      trackPiScore[iTrk] = prongCNN_out[0][0][3].item()
      trackPrScore[iTrk] = prongCNN_out[0][0][4].item()
      trackComp[iTrk] = prongCNN_out[1].item()
      trackPurity[iTrk] = prongCNN_out[2].item()
      trackProcess[iTrk] = prongCNN_out[3].argmax(1).item()
      trackPrimaryScore[iTrk] = prongCNN_out[3][0][0].item()
      trackFromNeutralScore[iTrk] = prongCNN_out[3][0][1].item()
      trackFromChargedScore[iTrk] = prongCNN_out[3][0][2].item()
      foundEnergy = True
      if trackMuScore[iTrk] >= trackPiScore[iTrk] and trackMuScore[iTrk] >= trackPrScore[iTrk]:
        trackRecoE[iTrk] = vertex.track_kemu_v[iTrk]
      elif trackPrScore[iTrk] >= trackMuScore[iTrk] and trackPrScore[iTrk] >= trackPiScore[iTrk]:
        trackRecoE[iTrk] = vertex.track_keproton_v[iTrk]
      else:
        try:
          trackRecoE[iTrk] = piKEestimator.Eval(getTrackLength(vertex.track_v[iTrk]))
        except:
          trackRecoE[iTrk] = -1.
          foundEnergy = False
      if foundEnergy:
        recoNuE[0] += trackRecoE[iTrk]

      if args.isMC:
        pdg, trackid, trueE, purity, completeness, allPdgs, allPurities = getMCProngParticle(prong_vv, mcpg, mcpm, adc_v, ioll)
        trackTruePID[iTrk] = pdg
        trackTrueTID[iTrk] = trackid
        trackTrueE[iTrk] = trueE
        trackTruePurity[iTrk] = purity
        trackTrueComp[iTrk] = completeness
        trackTrueElPurity[iTrk] = 0.
        trackTruePhPurity[iTrk] = 0.
        trackTrueMuPurity[iTrk] = 0.
        trackTruePiPurity[iTrk] = 0.
        trackTruePrPurity[iTrk] = 0.
        for iTru, current_pdg in enumerate(allPdgs):
          if current_pdg == 11:
            trackTrueElPurity[iTrk] = allPurities[iTru]
          if current_pdg == 22:
            trackTruePhPurity[iTrk] = allPurities[iTru]
          if current_pdg == 13:
            trackTrueMuPurity[iTrk] = allPurities[iTru]
          if current_pdg == 211:
            trackTruePiPurity[iTrk] = allPurities[iTru]
          if current_pdg == 2212:
            trackTruePrPurity[iTrk] = allPurities[iTru]
      #else:
      #  trackTruePID[iTrk] = 0
      #  trackTrueTID[iTrk] = -1
      #  trackTrueE[iTrk] = -1.
      #  trackTruePurity[iTrk] = -1.
      #  trackTrueComp[iTrk] = -1.
      #  trackTrueElPurity[iTrk] = -1.
      #  trackTruePhPurity[iTrk] = -1.
      #  trackTrueMuPurity[iTrk] = -1.
      #  trackTruePiPurity[iTrk] = -1.
      #  trackTruePrPurity[iTrk] = -1.
    #++++++ end track loop ++++++++++++++++++++++++++++++++++++++++++++++++++=


    #++++++ begin shower loop ++++++++++++++++++++++++++++++++++++++++++++++++++=
    for iShw, shower in enumerate(vertex.shower_v):

      for hit in shower:
        eventLarflowCluster.push_back(hit)
        if prongsAreContained:
          hitPt = rt.TVector3(hit[0], hit[1], hit[2])
          prongsAreContained = isFiducialWCSCE(hitPt)

      vertexNHits += shower.size()
      showerIsSecondary[iShw] = vertex.shower_isSecondary_v[iShw]
      showerNHits[iShw] = shower.size()
      showerCharge[iShw], vertexPixels, vertexCharge = addClusterCharge(iolcv,shower,vertexPixels,vertexCharge, 10.)
      showerCosTheta[iShw] = getCosThetaBeamShower(vertex.shower_trunk_v[iShw])
      showerCosThetaY[iShw] = getCosThetaGravShower(vertex.shower_trunk_v[iShw])
      showerDistToVtx[iShw] = getVertexDistance(vertex.shower_trunk_v[iShw].Vertex(), vertex)
      showerStartPosX[iShw] = vertex.shower_trunk_v[iShw].Vertex().X()
      showerStartPosY[iShw] = vertex.shower_trunk_v[iShw].Vertex().Y()
      showerStartPosZ[iShw] = vertex.shower_trunk_v[iShw].Vertex().Z()
      showerDir = getDirection(vertex.shower_trunk_v[iShw].Vertex(), vertex.shower_trunk_v[iShw].End())
      showerStartDirX[iShw], showerStartDirY[iShw], showerStartDirZ[iShw] = showerDir
      showerRecoE[iShw] = vertex.shower_plane_mom_vv[iShw][2].E()
      recoNuE[0] += vertex.shower_plane_mom_vv[iShw][2].E()

      cropPt = vertex.shower_trunk_v[iShw].Vertex()
      prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,shower,cropPt,10.,512,512)
      skip = False
      for p in range(3):
        if prong_vv[p].size() < 10:
          skip = True
          break
      if skip:
        showerClassified[iShw] = 0
        showerPID[iShw] = 0
        showerElScore[iShw] = -99.
        showerPhScore[iShw] = -99.
        showerMuScore[iShw] = -99.
        showerPiScore[iShw] = -99.
        showerPrScore[iShw] = -99.
        showerComp[iShw] = -1.
        showerPurity[iShw] = -1.
        showerProcess[iShw] = -1
        showerPrimaryScore[iShw] = -99.
        showerFromNeutralScore[iShw] = -99.
        showerFromChargedScore[iShw] = -99.
        continue

      with torch.no_grad():
        prongImage = makeImage(prong_vv).to(args.device)
        prongCNN_out = model(prongImage)
      showerClassified[iShw] = 1
      showerPID[iShw] = getPID(prongCNN_out[0].argmax(1).item())
      showerElScore[iShw] = prongCNN_out[0][0][0].item()
      showerPhScore[iShw] = prongCNN_out[0][0][1].item()
      showerMuScore[iShw] = prongCNN_out[0][0][2].item()
      showerPiScore[iShw] = prongCNN_out[0][0][3].item()
      showerPrScore[iShw] = prongCNN_out[0][0][4].item()
      showerComp[iShw] = prongCNN_out[1].item()
      showerPurity[iShw] = prongCNN_out[2].item()
      showerProcess[iShw] = prongCNN_out[3].argmax(1).item()
      showerPrimaryScore[iShw] = prongCNN_out[3][0][0].item()
      showerFromNeutralScore[iShw] = prongCNN_out[3][0][1].item()
      showerFromChargedScore[iShw] = prongCNN_out[3][0][2].item()

      if args.isMC:
        pdg, trackid, trueE, purity, completeness, allPdgs, allPurities = getMCProngParticle(prong_vv, mcpg, mcpm, adc_v, ioll)
        showerTruePID[iShw] = pdg
        showerTrueTID[iShw] = trackid
        showerTrueE[iShw] = trueE
        showerTruePurity[iShw] = purity
        showerTrueComp[iShw] = completeness
        showerTrueElPurity[iShw] = 0.
        showerTruePhPurity[iShw] = 0.
        showerTrueMuPurity[iShw] = 0.
        showerTruePiPurity[iShw] = 0.
        showerTruePrPurity[iShw] = 0.
        for iTru, current_pdg in enumerate(allPdgs):
          if current_pdg == 11:
            showerTrueElPurity[iShw] = allPurities[iTru]
          if current_pdg == 22:
            showerTruePhPurity[iShw] = allPurities[iTru]
          if current_pdg == 13:
            showerTrueMuPurity[iShw] = allPurities[iTru]
          if current_pdg == 211:
            showerTruePiPurity[iShw] = allPurities[iTru]
          if current_pdg == 2212:
            showerTruePrPurity[iShw] = allPurities[iTru]
      #else:
      #  showerTruePID[iShw] = 0
      #  showerTrueTID[iShw] = -1
      #  showerTrueE[iShw] = -1.
      #  showerTruePurity[iShw] = -1.
      #  showerTrueComp[iShw] = -1.
      #  showerTrueElPurity[iShw] = -1.
      #  showerTruePhPurity[iShw] = -1.
      #  showerTrueMuPurity[iShw] = -1.
      #  showerTruePiPurity[iShw] = -1.
      #  showerTruePrPurity[iShw] = -1.
    #++++++ end shower loop ++++++++++++++++++++++++++++++++++++++++++++++++++=

    if vtxIsFiducial[0] == 1:
      if prongsAreContained:
        vtxContainment[0] = 2
      else:
        vtxContainment[0] = 1

    for i in range(nTracks[0]):
      trackHitFrac[i] = trackNHits[i] / (1.0*vertexNHits)
      trackChargeFrac[i] = trackCharge[i] / vertexCharge

    for i in range(nShowers[0]):
      showerHitFrac[i] = showerNHits[i] / (1.0*vertexNHits)
      showerChargeFrac[i] = showerCharge[i] / vertexCharge

    #Do PCA for all hits attached to vertex
    try:
 
      eventCluster = larflow.reco.cluster_from_larflowcluster(eventLarflowCluster)
      for iPCA in range(3):
        eventPCAxis0[iPCA] = eventCluster.pca_axis_v[0][iPCA]
        eventPCAxis1[iPCA] = eventCluster.pca_axis_v[1][iPCA]
        eventPCAxis2[iPCA] = eventCluster.pca_axis_v[2][iPCA]
        eventPCEigenVals[iPCA] = eventCluster.pca_eigenvalues[iPCA]

      #project event 3D points to first PC axis with 0 at vertex projection
      vtxProj = [0.,0.,0.]
      ldist = 0.
      for c in range(3):
        ldist += (vertex.pos[c] - eventCluster.pca_center[c])*eventPCAxis0[c]
      for c in range(3):
        vtxProj[c] = eventCluster.pca_center[c] + ldist*eventPCAxis0[c]

      projDists = []
      avgPosPDistTimes = 0
      avgNegPDistTimes = 0
      nPosPDists = 0
      nNegPDists = 0
      for hit in eventLarflowCluster:
        projPt = [0.,0.,0.]
        ldist = 0.
        for c in range(3):
          ldist += (hit[c] - eventCluster.pca_center[c])*eventPCAxis0[c]
        if ldist > 0.:
          avgPosPDistTimes += hit.tick
          nPosPDists += 1
        else:
          avgNegPDistTimes += hit.tick
          nNegPDists += 1
        for c in range(3):
          projPt[c] = eventCluster.pca_center[c] + ldist*eventPCAxis0[c]
        projDist = 0.
        for c in range(3):
          projDist += (projPt[c] - vtxProj[c])**2
        projDists.append(sqrt(projDist))
      projDists.sort()

      #check if PCA axis is pointing in direction of increaasing hit times
      avgPosPDistTimes /= (1.0*nPosPDists)
      avgNegPDistTimes /= (1.0*nNegPDists)
      if avgPosPDistTimes > 0.:
        eventPCAxis0TSlope[0] = 1
      else:
        eventPCAxis0TSlope[0] = -1

      #calculate maximum point gap along PCA projection and maximum charge in between gaps
      maxGapFull = -1.
      maxGap90 = -1.
      maxGap80 = -1.
      maxGap70 = -1.
      maxGap60 = -1.
      maxCntD02 = -1.
      maxCntD04 = -1.
      maxCntD06 = -1.
      maxCntD08 = -1.
      maxCntD10 = -1.
      currentCntD02 = 0.
      currentCntD04 = 0.
      currentCntD06 = 0.
      currentCntD08 = 0.
      currentCntD10 = 0.

      for iEP in range(1, len(projDists)):

        gap = projDists[iEP] - projDists[iEP-1]

        if gap > maxGapFull:
          maxGapFull = gap
        if iEP < int(0.9*len(projDists)) and gap > maxGap90:
          maxGap90 = gap
        if iEP < int(0.8*len(projDists)) and gap > maxGap80:
          maxGap80 = gap
        if iEP < int(0.7*len(projDists)) and gap > maxGap70:
          maxGap70 = gap
        if iEP < int(0.6*len(projDists)) and gap > maxGap60:
          maxGap60 = gap

        if gap > 2. or iEP == (len(projDists) - 1):
          if currentCntD02 > maxCntD02:
            maxCntD02 = currentCntD02
          currentCntD02 = 0.
        else:
          currentCntD02 += gap
        if gap > 4. or iEP == (len(projDists) - 1):
          if currentCntD04 > maxCntD04:
            maxCntD04 = currentCntD04
          currentCntD04 = 0.
        else:
          currentCntD04 += gap
        if gap > 6. or iEP == (len(projDists) - 1):
          if currentCntD06 > maxCntD06:
            maxCntD06 = currentCntD06
          currentCntD06 = 0.
        else:
          currentCntD06 += gap
        if gap > 8. or iEP == (len(projDists) - 1):
          if currentCntD08 > maxCntD08:
            maxCntD08 = currentCntD08
          currentCntD08 = 0.
        else:
          currentCntD08 += gap
        if gap > 10. or iEP == (len(projDists) - 1):
          if currentCntD10 > maxCntD10:
            maxCntD10 = currentCntD10
          currentCntD10 = 0.
        else:
          currentCntD10 += gap

      eventPCProjMaxGap[0] = maxGapFull
      eventPCProjMaxGap[1] = maxGap90
      eventPCProjMaxGap[2] = maxGap80
      eventPCProjMaxGap[3] = maxGap70
      eventPCProjMaxGap[4] = maxGap60
      eventPCProjMaxDist[0] = maxCntD02
      eventPCProjMaxDist[1] = maxCntD04
      eventPCProjMaxDist[2] = maxCntD06
      eventPCProjMaxDist[3] = maxCntD08
      eventPCProjMaxDist[4] = maxCntD10

    except:
      print("warning: PCA failed")
      for iPCA in range(3):
        eventPCAxis0[iPCA] = 0.
        eventPCAxis1[iPCA] = 0.
        eventPCAxis2[iPCA] = 0.
        eventPCEigenVals[iPCA] = 0.
        eventPCAxis0TSlope[0] = 0
      for iPrj in range(5):
        eventPCProjMaxGap[iPrj] = -9.
        eventPCProjMaxDist[iPrj] = -9.

    eventTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


if args.isMC:
  totPOT[0] = totPOT_
  totGoodPOT[0] = totGoodPOT_
  potTree.Fill()

outRootFile.cd()
eventTree.Write("",rt.TObject.kOverwrite)
if args.isMC:
  potTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()

