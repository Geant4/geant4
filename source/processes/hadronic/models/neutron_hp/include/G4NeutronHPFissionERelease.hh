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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPFissionERelease.hh,v 1.6 2001-07-26 09:28:01 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFissionERelease_h
#define G4NeutronHPFissionERelease_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"

class G4NeutronHPFissionERelease
{
  public:
  G4NeutronHPFissionERelease(){}
  ~G4NeutronHPFissionERelease(){}
  
  inline void Init(G4std::ifstream & aDataFile)
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
            
    fragmentKinetic*=eV;
    promptNeutronKinetic*=eV;
    delayedNeutronKinetic*=eV;
    promptGammaEnergy*=eV;
    delayedGammaEnergy*=eV;
    delayedBetaEnergy*=eV;
    neutrinoEnergy*=eV;
    reducedTotalEnergy*=eV;
    totalEnergy*=eV;
  }
  
  inline G4double GetTotalEnergy(G4double deltaNNeu, G4double anEnergy) 
  {
     G4double result, delta, energy;
     energy = anEnergy/eV;
     delta = -(1.057*energy - 8.07*deltaNNeu);
     result = totalEnergy - delta*eV;
     return result;
  }
  inline G4double GetFragmentKinetic() 
  {
     return fragmentKinetic;
  }
  inline G4double GetPromptNeutronKinetic(G4double deltaNNeu, G4double anEnergy)
  {
     G4double result, delta, energy;
     energy = anEnergy/eV;
     delta = -(1.307*energy - 8.07*deltaNNeu);
     result = totalEnergy - delta*eV;
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
