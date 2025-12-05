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


#include "RunAction.hh"

#include "globals.hh"
#include "G4Threading.hh"

#include "G4IAEAphspWriter.hh"
#include "G4IAEAphspWriterStack.hh"
#include "IAEAphspRun.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  G4cout << "Destroying RunAction object" << G4endl;
  if (fIAEAphspWriterStack) {
    fIAEAphspWriterStack->ClearZphspVec();
    delete fIAEAphspWriterStack;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated to store run-persistent data.

  if ( G4Threading::IsMultithreadedApplication() ) {
    if (!(IsMaster()) && fIAEAphspWriterStack ) {
      G4cout << "Generating a worker IAEAphspRun with IAEAphspWriterStack!!"
	     << G4endl;
      return new IAEAphspRun(fIAEAphspWriterStack);
    }
    else {
      if (IsMaster()) G4cout << "Generating a master IAEAphspRun!!" << G4endl;
      else G4cout << "Generating a worker IAEAphspRun!!" << G4endl;
      return new IAEAphspRun();
    }
  }
  else {  // sequential mode
    if (fIAEAphspWriterStack)
      return new IAEAphspRun(fIAEAphspWriterStack);
    else
      return new IAEAphspRun();
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "RunAction::EndOfRunAction() " << G4endl;

  if ( G4Threading::IsMultithreadedApplication() ) {
    if (IsMaster()) {
      auto masterRun = static_cast<const IAEAphspRun*>(aRun);
      auto iaeaphspWriter = masterRun->GetIAEAphspWriter();
      if (iaeaphspWriter) {
	// The IAEAphsp files are open at first call of IAEAphspRun::Merge()
	iaeaphspWriter->CloseIAEAphspOutFiles();
      }
    }
  }
  else {    // sequential mode
    const IAEAphspRun* constRun = dynamic_cast<const IAEAphspRun*>(aRun);
    IAEAphspRun* iaeaRun = const_cast<IAEAphspRun*>(constRun);

    auto phspStack = iaeaRun->GetIAEAphspWriterStack();

    if (phspStack) {  // We defined IAEAphsp stack
      iaeaRun->DumpToIAEAphspFiles(phspStack);

      // Update the number of original histories to all files and close
      const G4int histories = iaeaRun->GetNumberOfEvent();
      const size_t nPhsp = phspStack->GetZphspVec()->size();
      auto iaeaphspWriter = iaeaRun->GetIAEAphspWriter();
      if (iaeaphspWriter) {
	for (size_t jj = 0; jj < nPhsp; jj++)
	  iaeaphspWriter->SumOrigHistories(jj, histories);

	iaeaphspWriter->CloseIAEAphspOutFiles();
      }
      else {
	G4ExceptionDescription msg;
	msg << "Could not get G4IAEAphspWriter object after dumping info"
	    << "from G4IAEAphspWriterStack."
	    << G4endl;
	G4Exception("RunAction::EndOfRunAction()",
		    "RunAction001", FatalException, msg);
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetIAEAphspWriterStack(const G4String& namePrefix)
{
  fIAEAphspWriterStack = new G4IAEAphspWriterStack(namePrefix);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddZphsp(const G4double val)
{
  if (fIAEAphspWriterStack) {
    fIAEAphspWriterStack->AddZphsp(val);
  }
  else {
    G4ExceptionDescription msg;
      msg << "z_phsp value passed, but no IAEAphsp output file name provided!"
	  << G4endl;
      G4Exception("RunAction::AddZphsp()",
		  "RunAction002", FatalException, msg);
  }
}
