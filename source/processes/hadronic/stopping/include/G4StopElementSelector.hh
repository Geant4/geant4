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
// $Id$
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 2 April 2000
//
// Class Description: 
//
// Selection of elements for negative particle capture
// Selection between decay/capture for mu-
// N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.

//-----------------------------------------------------------------------------
//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description
// 17/05/2006  V.Ivanchenko Cleanup
//
//-----------------------------------------------------------------------------

#ifndef G4StopElementSelector_h
#define G4StopElementSelector_h 1
 
#include "globals.hh"
#include "G4Element.hh"

class G4Material;

class G4StopElementSelector 
{ 
public:
 
  G4StopElementSelector();
  
  ~G4StopElementSelector();

  G4Element* GetElement(const G4Material* aMaterial);
  G4double   GetMuonCaptureRate(G4double Z, G4double A);
  G4double   GetMuonDecayRate(G4double Z, G4double A);

private:
  // hide assignment operator as private 
  G4StopElementSelector& operator=(const G4StopElementSelector &right);
  G4StopElementSelector(const G4StopElementSelector& );

};

#endif
 
