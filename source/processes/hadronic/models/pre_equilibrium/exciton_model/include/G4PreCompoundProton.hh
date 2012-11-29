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
// by V. Lara
//
// Modified:
// J. M. Quesada (July 08) cleanup 
// 20.08.2010 V.Ivanchenko added int Z and A and cleanup; added 
//                        G4ParticleDefinition to constructor,
//                        moved constructor and destructor to source

#ifndef G4PreCompoundProton_h
#define G4PreCompoundProton_h 1

#include "G4PreCompoundNucleon.hh"
#include "G4ProtonCoulombBarrier.hh"

class G4PreCompoundProton : public G4PreCompoundNucleon
{
public:

  G4PreCompoundProton();

  virtual ~G4PreCompoundProton();

protected:

  virtual G4double GetRj(G4int NumberParticles, G4int NumberCharged);

  virtual G4double CrossSection(G4double ekin) ; 

  virtual G4double GetAlpha();

  virtual G4double GetBeta();

  G4double GetOpt1(G4double K);

  G4double GetOpt2(G4double K);

  G4double GetOpt3(G4double K);

private:

  // operators
  G4PreCompoundProton(const G4PreCompoundProton &right);
  const G4PreCompoundProton& 
  operator= (const G4PreCompoundProton &right);
  G4int operator==(const G4PreCompoundProton &right) const;
  G4int operator!=(const G4PreCompoundProton &right) const;    

  G4ProtonCoulombBarrier theProtonCoulombBarrier;
  G4double ResidualAthrd;
  G4double FragmentAthrd;
  G4int FragmentA;
  G4int ResidualA;
  G4int ResidualZ;
  G4int theA;
  G4int theZ;
};
#endif
 
