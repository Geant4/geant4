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
///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hall.hh
// Description: Equipped to construct G4 geometry of the experimental hall
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Hall_h
#define CCalG4Hall_h 1

#include "CCalHall.hh"
#include "CCalG4Able.hh"

class CCalG4Hall: public CCalHall, public CCalG4Able {
public:
  //Constructor and Destructor
  CCalG4Hall(const G4String &name);
  virtual ~CCalG4Hall();
  
protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);
  virtual void constructDaughters();
};

#endif
