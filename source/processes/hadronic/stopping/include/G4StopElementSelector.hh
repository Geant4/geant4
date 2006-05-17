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
//
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
// **************************************************************
//
// File: G4StopElementSelector
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 2 April 2000
//
// Class Description: 
//
// Selection of elements for negative particle cupture
// Selection between decay/capture for mu-
// N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.
//
// Class Description: End 
//
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
  G4double  GetMuonCaptureRate(G4double Z, G4double A);
  G4double  GetMuonDecayRate(G4double Z, G4double A);

private:
  // hide assignment operator as private 
  G4StopElementSelector& operator=(const G4StopElementSelector &right);
  G4StopElementSelector(const G4StopElementSelector& );

};

#endif
 
