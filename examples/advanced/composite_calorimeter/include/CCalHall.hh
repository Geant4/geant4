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
// File: CCalHall.hh
// Description: Equipped to construct the geometry of the 96 Test Beam
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalHall_h
#define CCalHall_h 1

#include "CCalDetector.hh"

class CCalHall: public CCalDetector
{
public:
  //Constructor and Destructor
  CCalHall(const G4String &name);
  virtual ~CCalHall();

  //Get Methods
  G4String getMaterial()                    const {return genMaterial;}
  G4double   getDy_2Hall()                  const {return dy_2Hall;}
  G4double   getDx_2Hall()                  const {return dx_2Hall;}

protected:
  virtual G4int readFile();
  virtual void constructDaughters();

private:
  G4String genMaterial;            //General material
  G4double   dy_2Hall;               //Half width     of the Experimental Hall
  G4double   dx_2Hall;               //Half thickness of the Experimental Hall
};

#endif
