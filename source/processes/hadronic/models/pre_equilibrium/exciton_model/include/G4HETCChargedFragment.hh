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
// $Id: G4HETCChargedFragment.hh 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara

#ifndef G4HETCChargedFragment_h
#define G4HETCChargedFragment_h 1

#include "G4HETCFragment.hh"

class G4HETCChargedFragment : public G4HETCFragment
{
public:

  G4HETCChargedFragment(const G4ParticleDefinition*,
			G4VCoulombBarrier * aCoulombBarrier);

  virtual ~G4HETCChargedFragment();

  virtual G4double SampleKineticEnergy(const G4Fragment & aFragment);

private:

  // operators  
  G4HETCChargedFragment();
  G4HETCChargedFragment(const G4HETCChargedFragment &right);
  const G4HETCChargedFragment & 
  operator=(const G4HETCChargedFragment &right);
  G4bool operator==(const G4HETCChargedFragment &right) const; 
  G4bool operator!=(const G4HETCChargedFragment &right) const;    
    
};

#endif
