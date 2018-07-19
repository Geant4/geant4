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
#ifndef G4ParticleHPNBodyPhaseSpace_h
#define G4ParticleHPNBodyPhaseSpace_h 1

#include <fstream>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4Pow.hh"
#include "G4ios.hh"
#include "G4Neutron.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ReactionProduct.hh"

class G4ParticleHPNBodyPhaseSpace : public G4VParticleHPEnergyAngular
{
  public:
  
  G4ParticleHPNBodyPhaseSpace(){
   theTotalMass = 0.0;
   theTotalCount = 0;
  }
  ~G4ParticleHPNBodyPhaseSpace(){}
  
  public:
  
  void Init(G4double aMass, G4int aCount)
  {
    theTotalMass=aMass;
    theTotalCount=aCount;
  }

  void Init(std::istream & aDataFile)
  {
    aDataFile >> theTotalMass >> theTotalCount;
    theTotalMass *= G4Neutron::Neutron()->GetPDGMass();
  }
   
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
    
  private:
  
  inline G4double Prob(G4double anEnergy, G4double eMax, G4int n)
  {
    G4double result;
    result = std::sqrt(anEnergy)*G4Pow::GetInstance()->powA(eMax-anEnergy, 3.*n/2.-4.);
    return result;
  }
  
  inline G4double C(G4double anEnergy, G4double mass)
  {
    G4double result(0);
    if(theTotalCount==3) result = 4./CLHEP::pi/G4Pow::GetInstance()->powN(GetEmax(anEnergy, mass),2);
    if(theTotalCount==4) result = 105./32./G4Pow::GetInstance()->powA(GetEmax(anEnergy, mass), 3.5);
    //if(theTotalCount==5) result = 256./14./CLHEP::pi/G4Pow::GetInstance()->powA(GetEmax(anEnergy, mass), 5.);
    if(theTotalCount==5) result = 256./14./CLHEP::pi/G4Pow::GetInstance()->powN(GetEmax(anEnergy, mass), 5);
    return result;
  }
  
  inline G4double GetEmax(G4double anEnergy, G4double mass)
  {
    G4double result;
    G4double tMass = GetTarget()->GetMass();
    G4double pMass = GetProjectileRP()->GetMass();
    G4double availableEnergy = GetQValue() + anEnergy*tMass/(pMass+tMass);
    result = availableEnergy*(theTotalMass-mass)/theTotalMass;
    return result;
  }
  
  G4double MeanEnergyOfThisInteraction() {return -1; }
  
  private:
  
  G4double theTotalMass; 
  G4int theTotalCount;
  
};
#endif
