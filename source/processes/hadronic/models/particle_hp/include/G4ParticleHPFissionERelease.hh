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
// 070606 fix for Valgrind by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPFissionERelease_h
#define G4ParticleHPFissionERelease_h 1

#include <fstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"

class G4ParticleHPFissionERelease
{
  public:
     G4ParticleHPFissionERelease()
     : totalEnergy( 0.0 )
     , fragmentKinetic( 0.0 )
     , promptNeutronKinetic( 0.0 )
     , delayedNeutronKinetic( 0.0 )
     , promptGammaEnergy( 0.0 )
     , delayedGammaEnergy( 0.0 )
     , delayedBetaEnergy( 0.0 )
     , neutrinoEnergy( 0.0 )
     , reducedTotalEnergy( 0.0 )
     {
     }
  ~G4ParticleHPFissionERelease(){}
  
  inline void Init(std::istream & aDataFile)
  {
    G4double dummy;
    
    aDataFile >>dummy
              >>fragmentKinetic
              >>promptNeutronKinetic
              >>delayedNeutronKinetic
              >>promptGammaEnergy
              >>delayedGammaEnergy
              >>delayedBetaEnergy
              >>neutrinoEnergy
              >>reducedTotalEnergy
              >>totalEnergy;
            
    fragmentKinetic*=CLHEP::eV;
    promptNeutronKinetic*=CLHEP::eV;
    delayedNeutronKinetic*=CLHEP::eV;
    promptGammaEnergy*=CLHEP::eV;
    delayedGammaEnergy*=CLHEP::eV;
    delayedBetaEnergy*=CLHEP::eV;
    neutrinoEnergy*=CLHEP::eV;
    reducedTotalEnergy*=CLHEP::eV;
    totalEnergy*=CLHEP::eV;
  }
  
  inline G4double GetTotalEnergy(G4double deltaNNeu, G4double anEnergy) 
  {
     G4double result, delta, energy;
     energy = anEnergy/CLHEP::eV;
     delta = -(1.057*energy - 8.07*deltaNNeu);
     result = totalEnergy - delta*CLHEP::eV;
     return result;
  }
  inline G4double GetFragmentKinetic() 
  {
     return fragmentKinetic;
  }
  inline G4double GetPromptNeutronKinetic(G4double deltaNNeu, G4double anEnergy)
  {
     G4double result, delta, energy;
     energy = anEnergy/CLHEP::eV;
     delta = -(1.307*energy - 8.07*deltaNNeu);
     result = totalEnergy - delta*CLHEP::eV;
     return result;
  }
  inline G4double GetDelayedNeutronKinetic()
  {
    return delayedNeutronKinetic;
  }
  inline G4double GetPromptGammaEnergy()
  {
    return promptGammaEnergy;
  }
  inline G4double GetDelayedGammaEnergy(G4double anEnergy)
  {
    G4double delta = 0.075*anEnergy;
    G4double result = delayedGammaEnergy-delta;
    return result;
  }
  inline G4double GetDelayedBetaEnergy(G4double anEnergy)
  {
    G4double delta = 0.075*anEnergy;
    G4double result = delayedBetaEnergy-delta;
    return result;
  }
  inline G4double GetNeutrinoEnergy(G4double anEnergy)
  {
    G4double delta = 0.1*anEnergy;
    G4double result = neutrinoEnergy-delta;
    return result;
  }
  inline G4double GetReducedTotal(G4double deltaNNeu, G4double anEnergy)
  {
    return GetTotalEnergy(deltaNNeu, anEnergy) - GetNeutrinoEnergy(anEnergy);
  }
  private:
  
  G4double totalEnergy;
  G4double fragmentKinetic;
  G4double promptNeutronKinetic;
  G4double delayedNeutronKinetic;
  G4double promptGammaEnergy;
  G4double delayedGammaEnergy;
  G4double delayedBetaEnergy;
  G4double neutrinoEnergy;
  G4double reducedTotalEnergy;  // total - neutrino
};

#endif
