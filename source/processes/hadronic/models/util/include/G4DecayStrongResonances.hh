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
// $Id: G4DecayStrongResonances.hh 67984 2013-03-13 10:44:01Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4DecayStrongResonances
//
// Modified:  
// 02.11.2010 V.Ivanchenko moved constructor and destructor to source; 
//                         removed unused variable

#ifndef G4DecayStrongResonances_h
#define G4DecayStrongResonances_h 1

#include "G4KineticTrackVector.hh"
#include "G4ReactionProductVector.hh"

class G4V3DNucleus;

class G4DecayStrongResonances
{
public:
  G4DecayStrongResonances();
  ~G4DecayStrongResonances();

private:
   G4int operator==(G4DecayStrongResonances& right) {return (this == &right);}
   G4int operator!=(G4DecayStrongResonances& right) {return (this != &right);}
      
public:
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
                                      G4V3DNucleus* );
};

#endif // G4DecayStrongResonances_h


