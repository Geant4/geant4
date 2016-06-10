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
// $Id: G4PreCompoundAlpha.hh 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:
// J. M. Quesada (July 08) cleanup 
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor,
//                        moved constructor and destructor to source

#ifndef G4PreCompoundAlpha_h
#define G4PreCompoundAlpha_h 1

#include "G4PreCompoundIon.hh"
#include "G4AlphaCoulombBarrier.hh"

class G4PreCompoundAlpha : public G4PreCompoundIon
{
public:

  G4PreCompoundAlpha();

  virtual ~G4PreCompoundAlpha();

  virtual G4double GetRj(G4int NumberParticles, G4int NumberCharged) const;

  virtual G4double FactorialFactor(G4int N, G4int P) const;

  virtual G4double CoalescenceFactor(G4int A) const;

  virtual G4double GetAlpha() const;

private:

  // operators
  G4PreCompoundAlpha(const G4PreCompoundAlpha &right);
  const G4PreCompoundAlpha& 
  operator= (const G4PreCompoundAlpha &right);
  G4int operator==(const G4PreCompoundAlpha &right) const;
  G4int operator!=(const G4PreCompoundAlpha &right) const;    

  G4AlphaCoulombBarrier theAlphaCoulombBarrier;

};
#endif





 



