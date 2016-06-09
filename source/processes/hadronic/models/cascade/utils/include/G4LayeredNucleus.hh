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
//
// original class G4Nucleus by H.P. Wellisch
//   modified by J.L. Chuma, TRIUMF, 19-Nov-1996
//   last modified: 27-Mar-1997
//   Chr. Volcker, 10-Nov-1997: new methods and class variables.
//   M.G. Pia, 2 Oct 1998: modified GetFermiMomentum (original design was
//                       the source of memory leaks)
// This class G4LayeredNucleus by T. Lampen 14.6.2000

 
#ifndef G4LayeredNucleus_h
#define G4LayeredNucleus_h 1
 
#include "G4Nucleus.hh"
 
class G4LayeredNucleus  : public G4Nucleus
{
public:
  
  G4LayeredNucleus() : G4Nucleus()
  { 
  }
  
  G4LayeredNucleus( const G4double A, const G4double Z ) : G4Nucleus ( A,  Z )
  {
  }

  G4LayeredNucleus( const G4Material *aMaterial ) : G4Nucleus ( aMaterial )
  {
  }

  G4ThreeVector GetMomentum();
  void SetMomentum(const G4ThreeVector& mom);

private:
    
  G4ThreeVector momentumVector;

};
 
#endif
 
