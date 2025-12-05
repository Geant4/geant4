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


#include "IAEAphspRun.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4Run.hh"

#include <vector>

#include "G4IAEAphspWriter.hh"
#include "G4IAEAphspWriterStack.hh"


//==============================================================================

IAEAphspRun::IAEAphspRun()
:G4Run()
{
  G4cout << "Creating default IAEAphspRun object" << G4endl;
}



//==============================================================================

IAEAphspRun::IAEAphspRun(G4IAEAphspWriterStack* iaeaStack)
:G4Run()
{
  G4cout << "Creating IAEAphspRun object with IAEAphspWriterStack" << G4endl;
  fIAEAphspWriterStack = iaeaStack;
  fIAEAphspWriterStack->PrepareRun();
}



//==============================================================================

IAEAphspRun::~IAEAphspRun()
{
  G4cout << "Destroying IAEAphspRun object" << G4endl;

  if (fIAEAphspWriterStack)
    fIAEAphspWriterStack->ClearRunVectors();  // deletion only in RunAction!

  if (fIAEAphspWriter) delete fIAEAphspWriter;
}



//==============================================================================
//  RecordEvent() is a method called at end of event after EndOfEventAction().

void IAEAphspRun::RecordEvent(const G4Event* aEvent)
{
  G4Run::RecordEvent(aEvent);   // Mandatory to increment 'numberOfEvent'
  // G4cout << "IAEAphspRun: numberOfEvent = " << numberOfEvent << G4endl;
  // G4cout << "Event ID = " << aEvent->GetEventID() << G4endl;

  if (fIAEAphspWriterStack)
    fIAEAphspWriterStack->PrepareNextEvent();
}



//==============================================================================
// Merge info from local IAEAphspRun object to the global IAEAphspRun object

void IAEAphspRun::Merge(const G4Run* aRun)
{
  G4cout << "IAEAphspRun::Merge() started" << G4endl;

  const IAEAphspRun* localRun = static_cast<const IAEAphspRun*>(aRun);
  auto localPhspStack = localRun->GetIAEAphspWriterStack();

  if (localPhspStack) {     // only if we have IAEAphsp files
    DumpToIAEAphspFiles(localPhspStack);

    // Update the number of original histories to all files
    const G4int histories = localRun->GetNumberOfEvent();
    const size_t nPhsp = localPhspStack->GetZphspVec()->size();
    for (size_t jj = 0; jj < nPhsp; jj++)
      fIAEAphspWriter->SumOrigHistories(jj, histories);
  }

  G4Run::Merge(aRun);
}


//==============================================================================
// Dump info contained in local IAEAphspRun object into IAEAphsp files

void IAEAphspRun::DumpToIAEAphspFiles(const G4IAEAphspWriterStack* phspStack)
{
  // Get info from vectors in G4IAEAphspWriterStack
  if (phspStack) {
    auto localPdgMtrx   = phspStack->GetPDGMtrx();
    auto localPosMtrx   = phspStack->GetPosMtrx();
    auto localMomMtrx   = phspStack->GetMomMtrx();
    auto localEneMtrx   = phspStack->GetEneMtrx();
    auto localWtMtrx    = phspStack->GetWtMtrx();
    auto localNstatMtrx = phspStack->GetNstatMtrx();

    const size_t nPhsp = phspStack->GetZphspVec()->size();
    // -- DEBUG!!
    // G4cout << "IAEAphspRun: This run has " << nPhsp << " phsp planes stored."
    // 	   << G4endl;

    if (nPhsp != localPdgMtrx->size() || nPhsp != localPosMtrx->size() ||
	nPhsp != localMomMtrx->size() || nPhsp != localEneMtrx->size() ||
	nPhsp != localWtMtrx->size()  || nPhsp != localNstatMtrx->size() ) {
      G4ExceptionDescription msg;
      msg << "Number of zphsp stored != size of vectors storing phsp data."
	  << " MERGING IGNORED!" << G4endl;
      G4Exception("IAEAphspRun::DumpToIAEAphspFiles()",
		  "IAEAphspRun001", JustWarning, msg);
    }
    else {
      size_t jj = 0; // phsp plane counter
      for (const auto& phspPdgVec : *localPdgMtrx) {
	size_t nPart = phspPdgVec->size(); // Get number of particles
	// -- DEBUG!!
	// G4cout << "\tPhsp #" << jj << " stores " << nPart << " particles"
	//        << G4endl;

	if (nPart != (*localNstatMtrx)[jj]->size() ||
	    nPart != (*localPosMtrx)[jj]->size() ||
	    nPart != (*localMomMtrx)[jj]->size() ||
	    nPart != (*localEneMtrx)[jj]->size() ||
	    nPart != (*localWtMtrx)[jj]->size() ) {
	  G4ExceptionDescription msg;
	  msg << "Number of stored particles does not match in this "
	      << "thread-local run for phps plane #" << jj
	      << ". Merging ignored!" << G4endl;
	  G4Exception("IAEAphspRun::DumpToIAEAphspFiles()",
		      "IAEAphspRun002", JustWarning, msg);
	}
	else {
	  // Everything OK to dump particles into the IAEAphsp output files

	  // 1. If the G4IAEAphspWriter object was not created yet,
	  // create it, take data from G4IAEAphspWriterStack and open files
	  if (!fIAEAphspWriter) {
	    const G4String namePrefix = phspStack->GetFileName();
	    fIAEAphspWriter = new G4IAEAphspWriter(namePrefix);
	    fIAEAphspWriter->SetDataFromIAEAStack(phspStack);
	    fIAEAphspWriter->OpenIAEAphspOutFiles(this);
	  }

	  // 2. Loop over all vectors of this phsp to get the dynamic info
	  // and write it into the corresponding IAEAphsp file

	  size_t ii = 0;
	  for (const auto& pdg : *phspPdgVec) {
	    G4int nStat = (*(*localNstatMtrx)[jj])[ii];
	    G4double kinE = (*(*localEneMtrx)[jj])[ii];
	    G4double wt = (*(*localWtMtrx)[jj])[ii];
	    G4ThreeVector pos = (*(*localPosMtrx)[jj])[ii];
	    G4ThreeVector momDir = (*(*localMomMtrx)[jj])[ii];
	    fIAEAphspWriter->WriteIAEAParticle(jj, nStat, pdg, kinE, wt,
					       pos, momDir);
	    ii++;

	    // -- DEBUG!!
	    // G4cout << "\tPDG = " << pdg
	    // 	   << "  n_stat = " << (*(*localNstatMtrx)[jj])[ii]
	    // 	   << "  pos/cm = " << ((*(*localPosMtrx)[jj])[ii]) /cm
	    // 	   << "  momDir = " << (*(*localMomMtrx)[jj])[ii]
	    // 	   << "  kinE/MeV = " << ((*(*localEneMtrx)[jj])[ii]) /MeV
	    // 	   << "  wt = " << (*(*localWtMtrx)[jj])[ii] << G4endl;
	  }
	}
	jj++;
      }
    }
  }
  else {
    G4ExceptionDescription msg;
    msg << "This function is not meant to be called if no "
	<< "G4IAEAphspWriterStack object has been defined."
	<< G4endl;
    G4Exception("IAEAphspRun::DumpToIAEAphspFiles()",
		"IAEAphspRun003", FatalException, msg);
  }
}
