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
// File: G4RootEvent.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4RootEvent.hh"

// Addtional Include:
#include "G4Event.hh"
#include "CLHEP/HepMC/GenEvent.h"
#include "G4MCTEvent.hh"
#include "G4Pevent.hh"

ClassImp(G4RootEvent)

// Implementation of Destructor #1
G4RootEvent::~G4RootEvent()
{}

// Implementation of MakeTransientObject
G4Pevent* G4RootEvent::MakeTransientObject()
{
  HepMC::GenEvent* hepevt = 0;
  G4MCTEvent*        mctevt = 0;
  G4Event*         g4evt  = 0;
  // G4Pevent* anEvt = new G4Pevent(hepevt, mctevt, g4evt);
  G4Pevent* anEvt = new G4Pevent(hepevt, g4evt);
  return anEvt;
}

// End of G4RootEvent.cc

