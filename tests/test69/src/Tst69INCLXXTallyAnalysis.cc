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
#include "g4root.hh"
#include "Tst69INCLXXTallyAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadSecondary.hh"
#include "G4DynamicParticle.hh"
#include <sstream>

const char *Tst69INCLXXTallyAnalysis::incPartDict[nIncPart] = {
  "N",
  "Pi",
  "Alpha",
  "C12"
};

const char *Tst69INCLXXTallyAnalysis::outPartDict[nOutPart] = {
  "n",
  "p",
  "Pi",
  "Alpha",
  "Li7"
};

const char *Tst69INCLXXTallyAnalysis::incPartLongDict[nIncPart] = {
  "nucleon",
  "pion",
  "alpha",
  "C12"
};

const char *Tst69INCLXXTallyAnalysis::outPartLongDict[nOutPart] = {
  "neutrons",
  "protons",
  "pions",
  "alphas",
  "Li7"
};

const char *Tst69INCLXXTallyAnalysis::eSliceDict[nESlice] = {
  "LE",
  "ME",
  "HE"
};

const G4float Tst69INCLXXTallyAnalysis::eSliceLow[nESlice+1] = {
  0.,
  200.,
  3000.,
  20000
};

const G4int Tst69INCLXXTallyAnalysis::nBins = 60;

const G4int Tst69INCLXXTallyAnalysis::nZBins = 100;

Tst69INCLXXTallyAnalysis::Tst69INCLXXTallyAnalysis(G4String const &physList, const G4bool histos, const G4bool ntuple) :
  factoryOn(false),
  withHistos(histos),
  withNtuple(ntuple)
{
  filenameStem = "test69_tally_" + physList;
}

Tst69INCLXXTallyAnalysis::~Tst69INCLXXTallyAnalysis() {
  Close();
}

void Tst69INCLXXTallyAnalysis::Open() {
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  G4String extension = analysisManager->GetFileType();
  filename = filenameStem + "." + extension;

  // Open an output file
  G4bool fileOpen = analysisManager->OpenFile(filename);
  if (!fileOpen) {
    G4cout << "\n---> Tst69INCLXXTallyAnalysis::Open(): cannot open " << filename << G4endl;
    return;
  }

  // Create ntuple
  if(withNtuple) {
    analysisManager->SetFirstNtupleId(1);
    ntupleID = analysisManager->CreateNtuple("reac", "Reactions treated by INCL++");
    colIDAp = analysisManager->CreateNtupleIColumn("Ap");
    colIDZp = analysisManager->CreateNtupleIColumn("Zp");
    colIDEp = analysisManager->CreateNtupleFColumn("Ep");
    colIDAt = analysisManager->CreateNtupleIColumn("At");
    colIDZt = analysisManager->CreateNtupleIColumn("Zt");
    analysisManager->FinishNtuple();
  }

  // Create histos
  G4String hName, hDescription;
  std::stringstream ss;
  if(withHistos) {
    analysisManager->SetFirstHistoId(1);
    for(int iIncPart=0; iIncPart<nIncPart; ++iIncPart) {
      for(int iESlice=0; iESlice<nESlice; ++iESlice) {
        for(int iOutPart=0; iOutPart<nOutPart; ++iOutPart) {
          // create histograms for the ejectile energies
          hName = "h_";
          hName += incPartDict[iIncPart];
          hName += "_";
          hName += outPartDict[iOutPart];
          hName += "_";
          hName += eSliceDict[iESlice];
          ss.str("");
          ss << "Energy spectrum of "
            << outPartLongDict[iOutPart]
            << " from reactions induced by "
            << incPartLongDict[iIncPart]
            << " between "
            << eSliceLow[iESlice]
            << " and "
            << eSliceLow[iESlice+1]
            << " AMeV";
          hDescription = ss.str();

          G4cout << "Creating histogram " << hName << " with description: " << hDescription << G4endl;
          hID_ioe[iIncPart][iOutPart][iESlice] =
            analysisManager->CreateH1(hName, hDescription, nBins, 0., 1.1*eSliceLow[iESlice+1]);
        }

        // create histograms for the ejectile charges
        hName = "hZ_";
        hName += incPartDict[iIncPart];
        hName += "_";
        hName += eSliceDict[iESlice];
        ss.str("");
        ss << "Charge distribution of baryons produced in reactions induced by "
          << incPartLongDict[iIncPart]
            << " between "
            << eSliceLow[iESlice]
            << " and "
            << eSliceLow[iESlice+1]
            << " AMeV";
          hDescription = ss.str();

          G4cout << "Creating histogram " << hName << " with description: " << hDescription << G4endl;
          hID_Z_ie[iIncPart][iESlice] =
            analysisManager->CreateH1(hName, hDescription, nZBins, -0.5, nZBins-0.5);

      }
      // create histograms for the projectile energies
      hName = "h_";
      hName += incPartDict[iIncPart];
      ss.str("");
      ss << "Energy spectrum of reactions induced by "
        << incPartLongDict[iIncPart];
      hDescription = ss.str();

      G4cout << "Creating histogram " << hName << " with description: " << hDescription << G4endl;
      hID_i[iIncPart] =
        analysisManager->CreateH1(hName, hDescription, nBins, 0., 1.1*eSliceLow[nESlice]);

    }
  }

  factoryOn = true;
  G4cout << "\n----> Tree is opened in " << filename << G4endl;
}

void Tst69INCLXXTallyAnalysis::Close() {
  if(factoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    if(withHistos) {
      for(int iIncPart=0; iIncPart<nIncPart; ++iIncPart) {
        for(int iESlice=0; iESlice<nESlice; ++iESlice) {
          for(int iOutPart=0; iOutPart<nOutPart; ++iOutPart) {
            tools::histo::h1d *histo = analysisManager->GetH1(hID_ioe[iIncPart][iOutPart][iESlice]);
            const double mean = histo->mean();
            const double rms = histo->rms();
            G4cout
              << outPartLongDict[iOutPart] << " from "
              << incPartLongDict[iIncPart]
              << "-induced reactions (" << eSliceLow[iESlice]
              << "-" << eSliceLow[iESlice+1]
              << " MeV): " << mean << " +- " << rms << " MeV"
              << G4endl;
          }
        }
      }
    }

    analysisManager->Write();
    analysisManager->CloseFile();
    G4cout << "\n----> Tree is saved in " << filename << G4endl;

    delete G4AnalysisManager::Instance();
    factoryOn = false;
  }
}

void Tst69INCLXXTallyAnalysis::Tally(G4HadProjectile const &aTrack, G4Nucleus const &theNucleus, G4HadFinalState const &result) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4ParticleDefinition const * const trackDefinition = aTrack.GetDefinition();
  const G4int Ap = trackDefinition->GetAtomicMass();
  const G4int Zp = (G4int) trackDefinition->GetPDGCharge();
  const G4float Ep = aTrack.GetKineticEnergy()/MeV;
  const G4int At = theNucleus.GetA_asInt();
  const G4int Zt = theNucleus.GetZ_asInt();

  // Fill ntuple
  if(withNtuple) {
    analysisManager->FillNtupleIColumn(ntupleID, colIDAp, Ap);
    analysisManager->FillNtupleIColumn(ntupleID, colIDZp, Zp);
    analysisManager->FillNtupleFColumn(ntupleID, colIDEp, Ep);
    analysisManager->FillNtupleIColumn(ntupleID, colIDAt, At);
    analysisManager->FillNtupleIColumn(ntupleID, colIDZt, Zt);
    analysisManager->AddNtupleRow(ntupleID);
  }

  // Fill histos
  if(withHistos) {
    // Determine the incoming particle
    const G4int iIncPart = whichIncomingParticle(Ap, Zp);
    if(iIncPart<0)
      return;

    // Fill the histogram for the projectile energies
    const G4double EpPerNucleon = (Ap>0 ? Ep/Ap : Ep);
    analysisManager->FillH1(hID_i[iIncPart], EpPerNucleon);

    // Determine the energy slice
    const G4int iESlice = whichEnergySlice(EpPerNucleon);
    if(iESlice<0)
      return;

    // Loop over the outgoing particles
    const G4int nParticles = result.GetNumberOfSecondaries();
    for(G4int i=0; i<nParticles; ++i) {

      G4DynamicParticle const * const p = result.GetSecondary(i)->GetParticle();
      G4ParticleDefinition const * const pdef = p->GetDefinition();
      const G4int A = pdef->GetAtomicMass();
      const G4int Z = (G4int) pdef->GetPDGCharge();

      // Fill the histogram for the ejectile charges
      if(A>0)
        analysisManager->FillH1(hID_Z_ie[iIncPart][iESlice], Z);

      const G4int iOutPart = whichOutgoingParticle(A, Z);
      if(iOutPart<0)
        continue;

      const G4double EKin = p->GetKineticEnergy()/MeV;
      const G4double EKinPerNucleon = (A>0 ? EKin/A : EKin);
      // Fill the histogram for the outgoing particles
      analysisManager->FillH1(hID_ioe[iIncPart][iOutPart][iESlice], EKinPerNucleon);
    }
  }
}

G4int Tst69INCLXXTallyAnalysis::whichIncomingParticle(const G4int A, const G4int Z) {
  if(A==1)
    return 0;
  else if(A==0)
    return 1;
  else if(A==4 && Z==2)
    return 2;
  else if(A==12 && Z==6)
    return 3;
  else
    return -1;
}

G4int Tst69INCLXXTallyAnalysis::whichOutgoingParticle(const G4int A, const G4int Z) {
  if(A==1 && Z==0)
    return 0;
  else if(A==1 && Z==1)
    return 1;
  else if(A==0)
    return 2;
  else if(A==4 && Z==2)
    return 3;
  else if(A==7 && Z==3)
    return 4;
  else
    return -1;
}

G4int Tst69INCLXXTallyAnalysis::whichEnergySlice(const G4double E) {
  for(int iESlice=0; iESlice<nESlice; ++iESlice) {
    if(E>=eSliceLow[iESlice] && E<eSliceLow[iESlice+1])
      return iESlice;
  }
  return -1;
}

