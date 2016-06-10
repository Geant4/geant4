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
// $Id: G4ElementSelector.hh 66367 2012-12-18 09:18:08Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class header file 
//
// File name:  G4ElementSelector
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 21 April 2012 on base of G4StopElementSelector
//
// Class Description: 
//
// Selection of elements for slow negative particle capture
// Alternative selector should inherit from this class 
//
// N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//
//-----------------------------------------------------------------------------

#ifndef G4ElementSelector_h
#define G4ElementSelector_h 1
 
#include "globals.hh"
#include "G4Element.hh"
#include "G4Track.hh"
#include <vector>

class G4Nucleus;

class G4ElementSelector 
{ 
public:
 
  G4ElementSelector();
  
  virtual ~G4ElementSelector();

  virtual G4Element* SelectZandA(const G4Track& track, G4Nucleus*);

private:

  // hide assignment operator as private 
  G4ElementSelector& operator=(const G4ElementSelector &right);
  G4ElementSelector(const G4ElementSelector& );

  std::vector<G4double> prob;
};

#endif
 
