//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// File: G4HitRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4HitRootIO.hh"

// Addtional Include:
#include "G4VPHitsCollectionIO.hh"

// Implementation of Constructor #1
G4HitRootIO::G4HitRootIO()
{
  f_G4VPHitIO = (G4VPHitIO*) this;
}

// Implementation of Destructor #1
G4HitRootIO::~G4HitRootIO()
{}

// Implementation of GetHitRootIO
G4HitRootIO* G4HitRootIO::GetHitRootIO()
{
  G4HitRootIO* hio=0;
  if ( f_G4VPHitIO == 0 ) {
    hio = new G4HitRootIO;
  }
  return hio;
}

// Implementation of Store
bool G4HitRootIO::Store(const G4HCofThisEvent* hcevt)
{
  // cast for non-const transient collection of this event.
  G4HCofThisEvent* ahcevt = (G4HCofThisEvent*) hcevt;

  if ( m_verbose > 2 ) {
    std::cout << "G4HitRootIO: storing G4HCofThisEvent, # of collections: "
              << ahcevt->GetNumberOfCollections() << std::endl;
  }

  // Loop through hit collections
  for (int i=0; i< ahcevt->GetNumberOfCollections(); i++ ) {
    // Get hit collection for each sensitive detector
    G4VHitsCollection* hc = ahcevt->GetHC(i);
    if ( m_verbose > 2 ) {
      std::cout << "G4HitRootIO: collection#: " << i
                << ", HitsCollection: " << hc << "." << std::endl;
    }
    if (hc!=0) {
      // find out I/O manager for this hit collection
      G4VPHitsCollectionIO* hitIOman =
        (G4VPHitsCollectionIO*) f_catalog->GetHCIOmanager(hc->GetSDname());
      if ( hitIOman != 0) {
        if ( hitIOman->Store(hc) ) {
          if ( m_verbose > 2 ) {
            std::cout << "G4HitRootIO::Store() #" << i << " storing "
                      << hc->GetSDname() << "." << std::endl;
          }
        } else {
          std::cerr << "G4HitRootIO::Store() failed in storing "
                    << hc->GetSDname() << "." << std::endl;
        }
      } else {
        std::cerr << "G4HitRootIO::Store(): Hit Collection I/O manager is not"
                  << " registered for detector: '" << hc->GetSDname()
                  << "', collection name: '" << hc->GetName()
                  << "'." << std::endl;
        return false;
      }
    }
  }

  return true;
}

// Implementation of Retrieve
bool G4HitRootIO::Retrieve(G4HCofThisEvent*& hcevt)
{
  // new digits collections for this event
  hcevt = new G4HCofThisEvent;

  // Loop through the registered IO managers
  for (size_t i=0; i< f_catalog->NumberOfHCIOmanager(); i++ ) {

    // find out I/O manager for this hit collection
    G4VPHitsCollectionIO* hitIOman =
        (G4VPHitsCollectionIO*) f_catalog->GetHCIOmanager(i);
    if ( hitIOman != 0) {
      G4VHitsCollection* hc;
      if ( hitIOman->Retrieve(hc) ) {
        hcevt->AddHitsCollection(i, hc);
      }
    }
  }
  return true;
}

// End of G4HitRootIO.cc

