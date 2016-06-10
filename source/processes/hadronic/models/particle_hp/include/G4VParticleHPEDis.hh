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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4VParticleHPEDis_h
#define G4VParticleHPEDis_h 1

#include <fstream>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4ParticleHPVector.hh"

class G4VParticleHPEDis
{
  public:
  G4VParticleHPEDis()
  {
  }
  virtual ~G4VParticleHPEDis()
  {
  }
  
  virtual void Init(std::istream & theData) = 0;
  
  virtual G4double GetFractionalProbability(G4double anEnergy) = 0;
  
  virtual G4double Sample(G4double anEnergy) = 0; 
  
  private:
    
};

#endif
