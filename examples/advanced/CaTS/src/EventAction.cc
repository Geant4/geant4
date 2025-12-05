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
// CaTS (Calorimetry and Tracking Simulation)
//
// Authors: Hans Wenzel and Soon Yung Jun
//          (Fermi National Accelerator Laboratory)
//
// History: October 18th, 2021 : first implementation
// ********************************************************************
//
/// \file EventAction.cc
/// \brief Implementation of the CaTS::EventAction class

// Geant4 headers
#include "G4Event.hh"
#include "G4Types.hh"
#include "G4SDManager.hh"

// project headers:
#include "EventAction.hh"
#include "ConfigurationManager.hh"
#include "PhotonSD.hh"
#include "PhotonHit.hh"

#ifdef WITH_ROOT
#include "RootIO.hh"
#endif

#ifdef WITH_G4OPTICKS
#include "OpticksFlags.hh"
#include "G4Opticks.hh"
#include "G4OpticksHit.hh"
#include <string>
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* event)
{
  G4bool verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  if (verbose)
  {
    G4cout << "EventAction::EndOfEventAction Event:   " << event->GetEventID()
           << G4endl;
  }
#ifdef WITH_G4OPTICKS
  if (ConfigurationManager::getInstance()->isEnable_opticks())
  {
    G4Opticks* g4ok = G4Opticks::Get();
    g4ok->propagateOpticalPhotons(event->GetEventID());
    unsigned num_hits = g4ok->getNumHit();
    G4cout << "EndOfEventAction: num_hits: " << num_hits << G4endl;
    if (num_hits > 0)
    {
      G4HCtable* hctable = G4SDManager::GetSDMpointer()->GetHCtable();
      for (G4int i = 0; i < hctable->entries(); ++i)
      {
        G4String sdn = hctable->GetSDname(i);
        std::size_t found = sdn.find("PhotonDetector");
        if (found != std::string::npos)
	{
          PhotonSD* aSD =
            (PhotonSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdn);
          aSD->AddOpticksHits();
        }
      }
    }
    if (verbose)
    {
      G4cout << "**************************************************" << G4endl;
      G4cout << " EndOfEventAction: numphotons:   " << g4ok->getNumPhotons()
             << " Gensteps: " << g4ok->getNumGensteps()
             << "  Maxgensteps:  " << g4ok->getMaxGensteps() << G4endl;
      G4cout << " EndOfEventAction: num_hits: " << g4ok->getNumHit() << G4endl;
      G4cout << g4ok->dbgdesc() << G4endl;
    }
    g4ok->reset();
    if (verbose)
    {
      G4cout << "========================== After reset: " << G4endl;
      G4cout << " EndOfEventAction: numphotons:   " << g4ok->getNumPhotons()
             << " Gensteps: " << g4ok->getNumGensteps()
             << "  Maxgensteps:  " << g4ok->getMaxGensteps() << G4endl;
      G4cout << "EndOfEventAction: num_hits: " << g4ok->getNumHit() << G4endl;
      G4cout << g4ok->dbgdesc() << G4endl;
      G4cout << "**************************************************" << G4endl;
    }
  } // end isEnable_opticks
#endif  // end WITH_G4OPTICKS

#ifdef WITH_ROOT
  // Now we deal with the Geant4 Hit collections.
  if (ConfigurationManager::getInstance()->isWriteHits())
  {
    RootIO::GetInstance()->Write(event);
  }
#endif
}
