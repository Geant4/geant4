//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SBBremTable
//
// Author:        Mihaly Novak
//
// Creation date: 15.07.2018
//
// Modifications:
//
// -------------------------------------------------------------------
//
#include "G4SBBremTable.hh"

#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "Randomize.hh"

#include "G4String.hh"

#include "G4Log.hh"
#include "G4Exp.hh"

#include "zlib.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

G4SBBremTable::G4SBBremTable()
 : fMaxZet(-1), fNumElEnergy(-1), fNumKappa(-1), fUsedLowEenergy(-1.),
   fUsedHighEenergy(-1.), fLogMinElEnergy(-1.), fILDeltaElEnergy(-1.)
{}

G4SBBremTable::~G4SBBremTable()
{
  ClearSamplingTables();
}

void G4SBBremTable::Initialize(const G4double lowe, const G4double highe)
{
  fUsedLowEenergy  = lowe;
  fUsedHighEenergy = highe;
  BuildSamplingTables();
  InitSamplingTables();
//  Dump();
}

// run-time method that samples energy transferred to the emitted gamma photon
double G4SBBremTable::SampleEnergyTransfer(const G4double eekin,
                                           const G4double leekin,
                                           const G4double gcut,
                                           const G4double dielSupConst,
                                           const G4int    iZet,
                                           const G4int    matCutIndx,
                                           const G4bool   isElectron)
{
  static const G4double kAlpha2Pi = CLHEP::twopi*CLHEP::fine_structure_const;
  const G4double zet = (G4double)iZet;
  const G4int   izet = std::max(std::min(fMaxZet, iZet),1);
  G4double eGamma    = 0.;
  // this should be checked in the caller
  // if (eekin<=gcut) return kappa;
  const G4double lElEnergy     = leekin;
  const SamplingTablePerZ* stZ = fSBSamplingTables[izet];
  // get the gamma cut of this Z that corresponds to the current mat-cuts
  const std::size_t gamCutIndx = stZ->fMatCutIndxToGamCutIndx[matCutIndx];
  // gcut was not found: should never happen (only in verbose mode)
  if (gamCutIndx >= stZ->fNumGammaCuts || stZ->fGammaECuts[gamCutIndx]!=gcut) {
    G4String msg = " Gamma cut="+std::to_string(gcut) + " [MeV] was not found ";
    msg += "in case of Z = " + std::to_string(iZet) + ". ";
    G4Exception("G4SBBremTable::SampleEnergyTransfer()","em0X",FatalException,
                msg.c_str());
  }
  const G4double lGCut = stZ->fLogGammaECuts[gamCutIndx];
  // get the random engine
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  // find lower e- energy bin
  G4bool isCorner = false; // indicate that the lower edge e- energy < gam-gut
  G4bool isSimply = false; // simply sampling: isCorner+lower egde is selected
  G4int elEnergyIndx   = stZ->fMaxElEnergyIndx;
  // only if e- ekin is below the maximum value(use table at maximum otherwise)
  if (eekin < fElEnergyVect[elEnergyIndx]) {
    const G4double val = (lElEnergy-fLogMinElEnergy)*fILDeltaElEnergy;
    elEnergyIndx       = (G4int)val;
    G4double pIndxH    = val-elEnergyIndx;
    // check if we are at limiting case: lower edge e- energy < gam-gut
    if (fElEnergyVect[elEnergyIndx]<=gcut) {
      // recompute the probability of taking the higher e- energy bin()
      pIndxH   = (lElEnergy-lGCut)/(fLElEnergyVect[elEnergyIndx+1]-lGCut);
      isCorner = true;
    }
    //
    if (rndmEngine->flat()<pIndxH) {
      ++elEnergyIndx;      // take the table at the higher e- energy bin
    } else if (isCorner) { // take the table at the lower  e- energy bin
      // special sampling need to be done if lower edge e- energy < gam-gut:
      // actually, we "sample" from a table "built" at the gamm-cut i.e. delta
      isSimply = true;
    }
  }
  // should never happen under normal conditions but add protection
  if (!stZ->fTablesPerEnergy[elEnergyIndx]) {
    return 0.;
  }
  // Do the photon energy sampling:
  const STable *st =  stZ->fTablesPerEnergy[elEnergyIndx];
  const std::vector<G4double>& cVect = st->fCumCutValues;
  const std::vector<STPoint>&  pVect = st->fSTable;
  const G4double minVal = cVect[gamCutIndx];

  // should never happen under normal conditions but add protection
  if (minVal >= 1.) {
    return 0.;
  }
  // some transfomrmtion variables used in the looop
  const G4double lCurKappaC  = lGCut - leekin;
  const G4double lUsedKappaC = lGCut - fLElEnergyVect[elEnergyIndx];
  // dielectric (always) and e+ correction suppressions (if the primary is e+)
  G4double suppression = 1.;
  G4double rndm[2];
  // rejection loop starts here
  do {
    rndmEngine->flatArray(2, rndm);
    G4double kappa = 1.0;
    if (!isSimply) {
      const G4double cumRV = rndm[0]*(1.-minVal)+minVal;
      // find lower index of the values in the Cumulative Function: use linear
      // instead of binary search because it's faster in our case
      const G4int cumLIndx = LinSearch(pVect, fNumKappa, cumRV)-1;
//      const G4int cumLIndx = std::lower_bound( pVect.begin(), pVect.end(), cumRV,
//                                    [](const STPoint& p, const double& cumV) {
//                                    return p.fCum<cumV; } ) - pVect.begin() -1;
      const STPoint& stPL  = pVect[cumLIndx];
      const G4double pA    = stPL.fParA;
      const G4double pB    = stPL.fParB;
      const G4double cumL  = stPL.fCum;
      const G4double cumH  = pVect[cumLIndx+1].fCum;
      const G4double lKL   = fLKappaVect[cumLIndx];
      const G4double lKH   = fLKappaVect[cumLIndx+1];
      const G4double dm1   = (cumRV-cumL)/(cumH-cumL);
      const G4double dm2   = (1.+pA+pB)*dm1;
      const G4double dm3   = 1.+dm1*(pA+pB*dm1);
      // kappa sampled at E_i e- energy
      const G4double lKappa = lKL+dm2/dm3*(lKH-lKL);
      // transform lKappa to [log(gcut/eekin),0] form [log(gcut/E_i),0]
      kappa  = G4Exp(lKappa*lCurKappaC/lUsedKappaC);
    } else {
//      const G4double upLimit = std::min(1.*CLHEP::eV,eekin-gcut);
//      kappa = (gcut+rndm[0]*upLimit)/eekin;
      kappa = 1.-rndm[0]*(1.-gcut/eekin);
    }
    // compute the emitted photon energy: k
    eGamma = kappa*eekin;
    const G4double invEGamma = 1./eGamma;
    // compute dielectric suppression: 1/(1+[gk_p/k]^2)
    suppression = 1./(1.+dielSupConst*invEGamma*invEGamma);
    // add positron correction if particle is e+
    if (!isElectron) {
      const G4double e1     = eekin - gcut;
      const G4double iBeta1 =  (e1 + CLHEP::electron_mass_c2)
                              / std::sqrt(e1*(e1 + 2.*CLHEP::electron_mass_c2));
      const G4double e2     = eekin - eGamma;
      const G4double iBeta2 =  (e2 + CLHEP::electron_mass_c2)
                              / std::sqrt(e2*(e2 + 2.*CLHEP::electron_mass_c2));
      const G4double dum    = kAlpha2Pi*zet*(iBeta1 - iBeta2);
      suppression = (dum > -12.) ? suppression*G4Exp(dum) : 0.;
    }
  } while (rndm[1] > suppression);
  // end of rejection loop
  // return the sampled photon energy value k
  return eGamma;
}


void G4SBBremTable::BuildSamplingTables() {
  // claer
  ClearSamplingTables();
  LoadSTGrid();
  // First elements and gamma cuts data structures need to be built:
  // loop over all material-cuts and add gamma cut to the list of elements
  const G4ProductionCutsTable
  *thePCTable = G4ProductionCutsTable::GetProductionCutsTable();
  // a temporary vector to store one element
  std::vector<std::size_t> vtmp(1,0);
  std::size_t numMatCuts = thePCTable->GetTableSize();
  for (G4int imc=0; imc<(G4int)numMatCuts; ++imc) {
    const G4MaterialCutsCouple *matCut = thePCTable->GetMaterialCutsCouple(imc);
    if (!matCut->IsUsed()) {
      continue;
    }
    const G4Material*           mat = matCut->GetMaterial();
    const G4ElementVector* elemVect = mat->GetElementVector();
    const G4int              indxMC = matCut->GetIndex();
    const G4double gamCut = (*(thePCTable->GetEnergyCutsVector(0)))[indxMC];
    const std::size_t numElems = elemVect->size();
    for (std::size_t ielem=0; ielem<numElems; ++ielem) {
      const G4Element *elem = (*elemVect)[ielem];
      const G4int izet = std::max(std::min(fMaxZet, elem->GetZasInt()),1);
      if (!fSBSamplingTables[izet]) {
        // create data structure but do not load sampling tables yet: will be
        // loaded after we know what are the min/max e- energies where sampling
        // will be needed during the simulation for this Z
        // LoadSamplingTables(izet);
        fSBSamplingTables[izet] = new SamplingTablePerZ();
      }
      // add current gamma cut to the list of this element data (only if this
      // cut value is still not tehre)
      const std::vector<double> &cVect = fSBSamplingTables[izet]->fGammaECuts;
      std::size_t indx = std::find(cVect.cbegin(), cVect.cend(), gamCut)-cVect.cbegin();
      if (indx==cVect.size()) {
        vtmp[0] = imc;
        fSBSamplingTables[izet]->fGamCutIndxToMatCutIndx.push_back(vtmp);
        fSBSamplingTables[izet]->fGammaECuts.push_back(gamCut);
        fSBSamplingTables[izet]->fLogGammaECuts.push_back(G4Log(gamCut));
        ++fSBSamplingTables[izet]->fNumGammaCuts;
      } else {
        fSBSamplingTables[izet]->fGamCutIndxToMatCutIndx[indx].push_back(imc);
      }
    }
  }
}

void G4SBBremTable::InitSamplingTables() {
  const std::size_t numMatCuts = G4ProductionCutsTable::GetProductionCutsTable()
                            ->GetTableSize();
  for (G4int iz=1; iz<fMaxZet+1; ++iz) {
    SamplingTablePerZ* stZ = fSBSamplingTables[iz];
    if (!stZ) continue;
    // Load-in sampling table data:
    LoadSamplingTables(iz);
    // init data
    for (G4int iee=0; iee<fNumElEnergy; ++iee) {
      if (!stZ->fTablesPerEnergy[iee])
        continue;
      const G4double elEnergy = fElEnergyVect[iee];
      // 1 indicates that gamma production is not possible at this e- energy
      stZ->fTablesPerEnergy[iee]->fCumCutValues.resize(stZ->fNumGammaCuts,1.);
      // sort gamma cuts and other members accordingly
      for (std::size_t i=0; i<stZ->fNumGammaCuts-1; ++i) {
        for (std::size_t j=i+1; j<stZ->fNumGammaCuts; ++j) {
          if (stZ->fGammaECuts[j]<stZ->fGammaECuts[i]) {
            G4double dum0                   = stZ->fGammaECuts[i];
            G4double dum1                   = stZ->fLogGammaECuts[i];
            std::vector<std::size_t>   dumv = stZ->fGamCutIndxToMatCutIndx[i];
            stZ->fGammaECuts[i]             = stZ->fGammaECuts[j];
            stZ->fLogGammaECuts[i]          = stZ->fLogGammaECuts[j];
            stZ->fGamCutIndxToMatCutIndx[i] = stZ->fGamCutIndxToMatCutIndx[j];
            stZ->fGammaECuts[j]             = dum0;
            stZ->fLogGammaECuts[j]          = dum1;
            stZ->fGamCutIndxToMatCutIndx[j] = dumv;
          }
        }
      }
      // set couple indices to store the corresponding gamma cut index
      stZ->fMatCutIndxToGamCutIndx.resize(numMatCuts,-1);
      for (std::size_t i=0; i<stZ->fGamCutIndxToMatCutIndx.size(); ++i) {
        for (std::size_t j=0; j<stZ->fGamCutIndxToMatCutIndx[i].size(); ++j) {
          stZ->fMatCutIndxToGamCutIndx[stZ->fGamCutIndxToMatCutIndx[i][j]] = i;
        }
      }
      // clear temporary vector
      for (std::size_t i=0; i<stZ->fGamCutIndxToMatCutIndx.size(); ++i) {
        stZ->fGamCutIndxToMatCutIndx[i].clear();
      }
      stZ->fGamCutIndxToMatCutIndx.clear();
      //  init for each gamma cut that are below the e- energy
      for (std::size_t ic=0; ic<stZ->fNumGammaCuts; ++ic) {
        const G4double gamCut = stZ->fGammaECuts[ic];
        if (elEnergy>gamCut) {
          // find lower kappa index; compute the 'xi' i.e. cummulative value for
          // gamCut/elEnergy
          const G4double cutKappa = std::max(1.e-12, gamCut/elEnergy);
          const std::size_t iKLow = (cutKappa>1.e-12)
          ? std::lower_bound(fKappaVect.cbegin(), fKappaVect.cend(), cutKappa)
            - fKappaVect.cbegin() -1
          : 0;
          const STPoint* stpL = &(stZ->fTablesPerEnergy[iee]->fSTable[iKLow]);
          const STPoint* stpH = &(stZ->fTablesPerEnergy[iee]->fSTable[iKLow+1]);
          const G4double pA   = stpL->fParA;
          const G4double pB   = stpL->fParB;
          const G4double etaL = stpL->fCum;
          const G4double etaH = stpH->fCum;
          const G4double alph = G4Log(cutKappa/fKappaVect[iKLow])
                               /G4Log(fKappaVect[iKLow+1]/fKappaVect[iKLow]);
          const G4double dum  = pA*(alph-1.)-1.-pB;
          G4double val = etaL;
          if (alph==0.) {
            stZ->fTablesPerEnergy[iee]->fCumCutValues[ic] = val;
          } else {
            val = -(dum+std::sqrt(dum*dum-4.*pB*alph*alph))/(2.*pB*alph);
            val = val*(etaH-etaL)+etaL;
            stZ->fTablesPerEnergy[iee]->fCumCutValues[ic] = val;
          }
        }
      }
    }
  }
}

// should be called only from LoadSamplingTables(G4int) and once
void G4SBBremTable::LoadSTGrid() {
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path) {
    G4Exception("G4SBBremTable::LoadSTGrid()","em0006",
                FatalException, "Environment variable G4LEDATA not defined");
    return;
  }
  const G4String fname =  G4String(path) + "/brem_SB/SBTables/grid";
  std::ifstream infile(fname,std::ios::in);
  if (!infile.is_open()) {
    G4String msgc = "Cannot open file: " + fname;
    G4Exception("G4SBBremTable::LoadSTGrid()","em0006",
                FatalException, msgc.c_str());
    return;
  }
  // get max Z, # electron energies and # kappa values
  infile >> fMaxZet;
  infile >> fNumElEnergy;
  infile >> fNumKappa;
  // allocate space for the data and get them in:
  // (1.) first eletron energy grid
  fElEnergyVect.resize(fNumElEnergy);
  fLElEnergyVect.resize(fNumElEnergy);
  for (G4int iee=0; iee<fNumElEnergy; ++iee) {
    G4double  dum;
    infile >> dum;
    fElEnergyVect[iee]  = dum*CLHEP::MeV;
    fLElEnergyVect[iee] = G4Log(fElEnergyVect[iee]);
  }
  // (2.) then the kappa grid
  fKappaVect.resize(fNumKappa);
  fLKappaVect.resize(fNumKappa);
  for (G4int ik=0; ik<fNumKappa; ++ik) {
    infile >> fKappaVect[ik];
    fLKappaVect[ik] = G4Log(fKappaVect[ik]);
  }
  // (3.) set size of the main container for sampling tables
  fSBSamplingTables.resize(fMaxZet+1,nullptr);
  // init electron energy grid related variables: use accurate values !!!
//  fLogMinElEnergy   = G4Log(fElEnergyVect[0]);
//  fILDeltaElEnergy  = 1./G4Log(fElEnergyVect[1]/fElEnergyVect[0]);
  const G4double elEmin  = 100.0*CLHEP::eV; //fElEnergyVect[0];
  const G4double elEmax  =  10.0*CLHEP::GeV;//fElEnergyVect[fNumElEnergy-1];
  fLogMinElEnergy  = G4Log(elEmin);
  fILDeltaElEnergy = 1./(G4Log(elEmax/elEmin)/(fNumElEnergy-1.0));
  // reset min/max energies if needed
  fUsedLowEenergy  = std::max(fUsedLowEenergy ,elEmin);
  fUsedHighEenergy = std::min(fUsedHighEenergy,elEmax);
  //
  infile.close();
}

void G4SBBremTable::LoadSamplingTables(G4int iz) {
  // check if grid needs to be loaded first
  if (fMaxZet<0) {
    LoadSTGrid();
  }
  // load data for a given Z only once
  iz = std::max(std::min(fMaxZet, iz),1);
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path) {
    G4Exception("G4SBBremTable::LoadSamplingTables()","em0006",
                FatalException, "Environment variable G4LEDATA not defined");
    return;
  }
  const G4String fname =  G4String(path) + "/brem_SB/SBTables/sTableSB_"
                        + std::to_string(iz);
  std::istringstream infile(std::ios::in);
  // read the compressed data file into the stream
  ReadCompressedFile(fname, infile);
  // the SamplingTablePerZ object was already created, set size of containers
  // then load sampling table data for each electron energies
  SamplingTablePerZ* zTable = fSBSamplingTables[iz];
  //
  // Determine min/max elektron kinetic energies and indices
  const G4double minGammaCut = zTable->fGammaECuts[ std::min_element(
                 std::cbegin(zTable->fGammaECuts),std::cend(zTable->fGammaECuts))
                -std::cbegin(zTable->fGammaECuts)];
  const G4double elEmin = std::max(fUsedLowEenergy, minGammaCut);
  const G4double elEmax = fUsedHighEenergy;
  // find low/high elecrton energy indices where tables will be needed
  // low:
  zTable->fMinElEnergyIndx = 0;
  if (elEmin>=fElEnergyVect[fNumElEnergy-1]) {
    zTable->fMinElEnergyIndx = fNumElEnergy-1;
  } else {
    zTable->fMinElEnergyIndx = G4int(std::lower_bound(fElEnergyVect.cbegin(),
                                                fElEnergyVect.cend(), elEmin)
                                     - fElEnergyVect.cbegin() -1);
  }
  // high:
  zTable->fMaxElEnergyIndx = 0;
  if (elEmax>=fElEnergyVect[fNumElEnergy-1]) {
    zTable->fMaxElEnergyIndx = fNumElEnergy-1;
  } else {
    // lower + 1
    zTable->fMaxElEnergyIndx = G4int(std::lower_bound(fElEnergyVect.cbegin(),
                                                fElEnergyVect.cend(), elEmax)
                                     - fElEnergyVect.cbegin());
  }
  // protect
  if (zTable->fMaxElEnergyIndx<=zTable->fMinElEnergyIndx) {
    return;
  }
  // load sampling tables that are needed: file is already in the stream
  zTable->fTablesPerEnergy.resize(fNumElEnergy, nullptr);
  for (G4int iee=0; iee<fNumElEnergy; ++iee) {
    // go over data that are not needed
    if (iee<zTable->fMinElEnergyIndx || iee>zTable->fMaxElEnergyIndx) {
      for (G4int ik=0; ik<fNumKappa; ++ik) {
        G4double dum;
        infile >> dum; infile >> dum; infile >> dum;
      }
    } else { // load data that are needed
      zTable->fTablesPerEnergy[iee] = new STable();
      zTable->fTablesPerEnergy[iee]->fSTable.resize(fNumKappa);
      for (G4int ik=0; ik<fNumKappa; ++ik) {
        STPoint &stP = zTable->fTablesPerEnergy[iee]->fSTable[ik];
        infile >> stP.fCum;
        infile >> stP.fParA;
        infile >> stP.fParB;
      }
    }
  }
}

// clean away all sampling tables and make ready to re-install
void G4SBBremTable::ClearSamplingTables() {
  for (G4int iz=0; iz<fMaxZet+1; ++iz) {
    if (fSBSamplingTables[iz]) {
      for (G4int iee=0; iee<fNumElEnergy; ++iee) {
        if (fSBSamplingTables[iz]->fTablesPerEnergy[iee]) {
          fSBSamplingTables[iz]->fTablesPerEnergy[iee]->fSTable.clear();
          fSBSamplingTables[iz]->fTablesPerEnergy[iee]->fCumCutValues.clear();
        }
      }
      fSBSamplingTables[iz]->fTablesPerEnergy.clear();
      fSBSamplingTables[iz]->fGammaECuts.clear();
      fSBSamplingTables[iz]->fLogGammaECuts.clear();
      fSBSamplingTables[iz]->fMatCutIndxToGamCutIndx.clear();
      //
      delete fSBSamplingTables[iz];
      fSBSamplingTables[iz] = nullptr;
    }
  }
  fSBSamplingTables.clear();
  fElEnergyVect.clear();
  fLElEnergyVect.clear();
  fKappaVect.clear();
  fLKappaVect.clear();
  fMaxZet = -1;
}

//void G4SBBremTable::Dump() {
//  G4cerr<< "\n  =====   Dumping ===== \n" << G4endl;
//  for (G4int iz=0; iz<fMaxZet+1; ++iz) {
//    if (fSBSamplingTables[iz]) {
//      G4cerr<< "   ----> There are " << fSBSamplingTables[iz]->fNumGammaCuts
//            << " g-cut for Z = " << iz << G4endl;
//      for (std::size_t ic=0; ic<fSBSamplingTables[iz]->fGammaECuts.size(); ++ic)
//        G4cerr<< "        i = " << ic << "  "
//              << fSBSamplingTables[iz]->fGammaECuts[ic] << G4endl;
//    }
//  }
//}

// find lower bin index of value: used in acse of CDF values i.e. val in [0,1)
// while vector elements in [0,1]
G4int G4SBBremTable::LinSearch(const std::vector<STPoint>& vect,
                               const G4int size,
                               const G4double val) {
  G4int i= 0;
  while (i + 3 < size) {
    if (vect [i + 0].fCum > val) return i + 0;
    if (vect [i + 1].fCum > val) return i + 1;
    if (vect [i + 2].fCum > val) return i + 2;
    if (vect [i + 3].fCum > val) return i + 3;
    i += 4;
  }
  while (i < size) {
    if (vect [i].fCum > val)
      break;
    ++i;
  }
  return i;
}

// uncompress one data file into the input string stream
void G4SBBremTable::ReadCompressedFile(const G4String &fname,
                                       std::istringstream &iss) {
  std::string *dataString = nullptr;
  std::string compfilename(fname+".z");
  // create input stream with binary mode operation and positioning at the end
  // of the file
  std::ifstream in(compfilename, std::ios::binary | std::ios::ate);
  if (in.good()) {
     // get current position in the stream (was set to the end)
     std::streamoff fileSize = in.tellg();
     // set current position being the beginning of the stream
     in.seekg(0,std::ios::beg);
     // create (zlib) byte buffer for the data
     Bytef *compdata = new Bytef[fileSize];
     while(in) {
        in.read((char*)compdata, fileSize);
     }
     // create (zlib) byte buffer for the uncompressed data
     uLongf complen    = (uLongf)(fileSize*4);
     Bytef *uncompdata = new Bytef[complen];
     while (Z_OK!=uncompress(uncompdata, &complen, compdata, fileSize)) {
        // increase uncompressed byte buffer
        delete[] uncompdata;
        complen   *= 2;
        uncompdata = new Bytef[complen];
     }
     // delete the compressed data buffer
     delete [] compdata;
     // create a string from the uncompressed data (will be deleted by the caller)
     dataString = new std::string((char*)uncompdata, (long)complen);
     // delete the uncompressed data buffer
     delete [] uncompdata;
  } else {
    std::string msg = "  Problem while trying to read "
                      + compfilename + " data file.\n";
    G4Exception("G4SBBremTable::ReadCompressedFile","em0006",
                FatalException,msg.c_str());
    return;
  }
  // create the input string stream from the data string
  if (dataString) {
    iss.str(*dataString);
    in.close();
    delete dataString;
  }
}
