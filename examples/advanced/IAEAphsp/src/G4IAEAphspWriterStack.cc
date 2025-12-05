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
// Author: M.A. Cortes-Giraldo
//
// 2025-08-27: Its objects work in local runs, their mission is to store the
//   info of the particles to be written in the IAEAphsp output files at
//   the end of the run. In other words, they constitute a stack for
//   IAEAphsp particles until the run finishes, when the IAEAphsp file is
//   actually written.
//


#include "G4IAEAphspWriterStack.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4IAEAphspWriter.hh"

#include <set>
#include <vector>


//==============================================================================

G4IAEAphspWriterStack::G4IAEAphspWriterStack(const G4String filename)
{
  fFileName = filename;
  fZphspVec = new std::vector<G4double>;
  fIncrNumberVec = new std::vector<G4int>;
  fPassingTracksVec = new std::vector< std::set<G4int>* >;
  fPDGMtrx = new std::vector< std::vector<G4int>* >;
  fPosMtrx = new std::vector< std::vector<G4ThreeVector>* >;
  fMomMtrx = new std::vector< std::vector<G4ThreeVector>* >;
  fEneMtrx = new std::vector< std::vector<G4double>* >;
  fWtMtrx = new std::vector< std::vector<G4double>* >;
  fNstatMtrx = new std::vector< std::vector<G4int>* >;

  G4cout << "G4IAEAphspWriterStack object constructed for files \""
	 << fFileName << "_<zphsp>cm.IAEA*\"" << G4endl;
}


//==============================================================================

G4IAEAphspWriterStack::~G4IAEAphspWriterStack()
{
  if (fZphspVec) delete fZphspVec;
  if (fIncrNumberVec)    delete fIncrNumberVec;
  if (fPassingTracksVec) delete fPassingTracksVec;
  if (fPDGMtrx)    delete fPDGMtrx;
  if (fPosMtrx)    delete fPosMtrx;
  if (fMomMtrx)    delete fMomMtrx;
  if (fEneMtrx)    delete fEneMtrx;
  if (fWtMtrx)     delete fWtMtrx;
  if (fNstatMtrx)  delete fNstatMtrx;
}


//==============================================================================

void G4IAEAphspWriterStack::AddZphsp(const G4double zphsp)
{
  fZphspVec->push_back(zphsp);
  G4cout << "G4IAEAphspWriterStack: Registered phase-space plane at z = "
	 << zphsp/cm << " cm." << G4endl;
}


//==============================================================================

void G4IAEAphspWriterStack::ClearZphspVec()
{
  G4cout << "G4IAEAphspWriterStack: Removing all registered phase-space planes!"
	 << G4endl;
  fZphspVec->clear();
}


//==============================================================================

void G4IAEAphspWriterStack::SetDataFromWriter(const G4IAEAphspWriter* writer)
{
  if (!writer) {
    G4ExceptionDescription msg;
    msg <<  "No G4IAEAphspWriter has been constructed!" << G4endl;
    G4Exception("G4IAEAphspWriterStack::SetDataFromWriter()",
		"IAEAphspWriterStack001", FatalException, msg );
    return;
  }
  else {
    fFileName = writer->GetFileName();

    if (writer->GetZphspVec()->size() > 0) {
      (*fZphspVec) = *(writer->GetZphspVec()); // copy objects, not pointers

      G4cout << "G4IAEAphspWriterStack::fFileName = " << fFileName << G4endl;
      G4cout << "G4IAEAphspWriterStack::fZphspVec->size() = "
	     << fZphspVec->size() << G4endl;
    }
    else {
      G4ExceptionDescription msg;
      msg << "No phsp plane z-coordinate has been defined!" << G4endl;
      G4Exception("G4IAEAphspWriterStack::SetDataFromWriter()",
		  "IAEAphspWriterStack002", FatalErrorInArgument, msg );
      return;
    }
  }
}


//==============================================================================

void G4IAEAphspWriterStack::PrepareRun()
{
  size_t nZphsps = fZphspVec->size();

  fIncrNumberVec->reserve(nZphsps);
  fPassingTracksVec->reserve(nZphsps);
  fPDGMtrx->reserve(nZphsps);
  fPosMtrx->reserve(nZphsps);
  fMomMtrx->reserve(nZphsps);
  fEneMtrx->reserve(nZphsps);
  fWtMtrx->reserve(nZphsps);
  fNstatMtrx->reserve(nZphsps);
  
  for (size_t ii = 0; ii < nZphsps; ii++) {
    fIncrNumberVec->push_back(0);

    auto aSet = new std::set<G4int>;
    fPassingTracksVec->push_back(aSet);
    auto pdgVec = new std::vector<G4int>;
    fPDGMtrx->push_back(pdgVec);
    auto posVec = new std::vector<G4ThreeVector>;
    fPosMtrx->push_back(posVec);
    auto momVec = new std::vector<G4ThreeVector>;
    fMomMtrx->push_back(momVec);
    auto eneVec = new std::vector<G4double>;
    fEneMtrx->push_back(eneVec);
    auto wtVec = new std::vector<G4double>;
    fWtMtrx->push_back(wtVec);
    auto nstatVec = new std::vector<G4int>;
    fNstatMtrx->push_back(nstatVec);
  }
  G4cout << "G4IAEAphspWriterStack::PrepareRun() done!" << G4endl;
}


//==============================================================================

void G4IAEAphspWriterStack::PrepareNextEvent()
{
  // Update all the incremental history numbers.
  for ( auto& ii : (*fIncrNumberVec) )
    ii++;

  // Remove the track ID's stored during this event.
  for ( auto& trackIDs : (*fPassingTracksVec) )
    trackIDs->clear();

  // -- DEBUG!!
  // G4cout << "G4IAEAphspWriterStack ready for the next event!" << G4endl;
}



//==============================================================================

void G4IAEAphspWriterStack::StoreParticleIfEligible(const G4Step* aStep)
{
  const G4ThreeVector postR = aStep->GetPostStepPoint()->GetPosition();
  const G4ThreeVector preR = aStep->GetPreStepPoint()->GetPosition();
  const G4double postZ = postR.z();
  const G4double preZ = preR.z();

  // Check what phsp planes are being crossed
  size_t phspIdx = 0;
  for (const auto& phspZ : (*fZphspVec) ) {
    if ( (postZ-phspZ)*(preZ-phspZ) < 0 ) {
      // Get trackID and check if it is already in fPassingTracksVec[ii]
      const G4int trackID = aStep->GetTrack()->GetTrackID();
      std::set<G4int>::iterator is;
      is = (*fPassingTracksVec)[phspIdx]->find(trackID);

      if ( is == (*fPassingTracksVec)[phspIdx]->end() ) {
	// This particle has not crossed this phsp plane before.
	// Then, put it on the stack if it is of a type foreseen
	// by the IAEAphsp format
	const G4int pdgCode =
	  aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
	if (pdgCode == 22 || pdgCode == 11 || pdgCode == -11 ||
	    pdgCode == 2112 || pdgCode == 2212)
	  StoreIAEAParticle(aStep, phspIdx, pdgCode);
      }
    }
    phspIdx++;
  }
}


//==============================================================================

void G4IAEAphspWriterStack::StoreIAEAParticle(const G4Step* aStep,
					      const G4int phspIndex,
					      const G4int pdgCode)
{
  const G4double zStop = (*fZphspVec)[phspIndex];
  const G4Track* aTrack = aStep->GetTrack();

  // Get step info
  // --------------------------------
  const G4ThreeVector postR = aStep->GetPostStepPoint()->GetPosition();
  const G4ThreeVector preR = aStep->GetPreStepPoint()->GetPosition();
  const G4double postZ = postR.z();
  const G4double preZ = preR.z();

  // Set kinetic energy
  G4double kinEnergy;
  if (pdgCode == 22 || pdgCode == 2112) {  // gamma or neutron
    kinEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  }
  else if (pdgCode == 11 || pdgCode == -11 ||
	   pdgCode == 2212 ) {  // electron, positron or proton
    const G4double postE = aStep->GetPostStepPoint()->GetKineticEnergy();
    const G4double preE = aStep->GetPreStepPoint()->GetKineticEnergy();
    kinEnergy = preE + (postE-preE)*(zStop-preZ)/(postZ-preZ);
  }
  else { // not a particle for the IAEA format
    G4ExceptionDescription ED;
    ED << "\"" << aTrack->GetDefinition()->GetParticleName()
       << "\" is not supported by the IAEAphsp format; not recorded.";
    G4Exception("G4IAEAphspWriterStack::StoreIAEAParticle()",
		"IAEAphspWriterStack003", JustWarning, ED );
    return;
  }

  // Position
  const G4ThreeVector phspPos = preR + (postR-preR)*(zStop-preZ)/(postZ-preZ);

  // Momentum direction
  const G4ThreeVector phspMomDir =
    aStep->GetPreStepPoint()->GetMomentumDirection();

  // Track weight
  const G4double wt = aTrack->GetWeight();

  // n_stat value
  const G4int nStat = (*fIncrNumberVec)[phspIndex];

  // Store info in stacking vectors
  // ------------------------------
  ((*fPDGMtrx)[phspIndex])->push_back(pdgCode);
  ((*fNstatMtrx)[phspIndex])->push_back(nStat);
  ((*fPosMtrx)[phspIndex])->push_back(phspPos);
  ((*fMomMtrx)[phspIndex])->push_back(phspMomDir);
  ((*fEneMtrx)[phspIndex])->push_back(kinEnergy);
  ((*fWtMtrx)[phspIndex])->push_back(wt);

  // Once stored, reset the incremental history number (n_stat = 0)
  (*fIncrNumberVec)[phspIndex] = 0;

  // And now register this trackID to protect against multiple crossers
  (*fPassingTracksVec)[phspIndex]->insert( aTrack->GetTrackID() );

  // -- DEBUG!!
  // G4cout << "G4IAEAphspWriterStack: Particle stored in phsp plane ["
  // 	 << phspIndex << "] at place #" << ((*fPDGMtrx)[phspIndex])->size()
  // 	 << " with the following values:" << G4endl;
  // G4cout << "\tPDG = " << pdgCode << "   nStat = " << nStat
  // 	 << "   phspPos = " << phspPos << "   phspMomDir = " << phspMomDir
  // 	 << "   kinEnergy = " << kinEnergy << "   weight = " << wt
  // 	 << G4endl;

}


//==============================================================================

void G4IAEAphspWriterStack::ClearRunVectors()
{
  // Clear the run vectors at run termination

  fIncrNumberVec->clear();

  for (auto& trackIDs : (*fPassingTracksVec) ) delete trackIDs;
  fPassingTracksVec->clear();

  for (auto& pdgVec : (*fPDGMtrx) ) delete pdgVec;
  fPDGMtrx->clear();

  for (auto& posVec : (*fPosMtrx) ) delete posVec;
  fPosMtrx->clear();

  for (auto& momVec : (*fMomMtrx) ) delete momVec;
  fMomMtrx->clear();

  for (auto& eneVec : (*fEneMtrx) ) delete eneVec;
  fEneMtrx->clear();

  for (auto& wtVec : (*fWtMtrx) ) delete wtVec;
  fWtMtrx->clear();

  for (auto& nStatVec : (*fNstatMtrx) ) delete nStatVec;
  fNstatMtrx->clear();

  G4cout << "G4IAEAphspWriterStack run vectors cleaned!" << G4endl;
}
