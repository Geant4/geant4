// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPNBodyPhaseSpace.hh,v 1.1 1999-01-07 16:13:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPNBodyPhaseSpace_h
#define G4NeutronHPNBodyPhaseSpace_h 1

#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"
#include "G4Neutron.hh"
#include "G4VNeutronHPEnergyAngular.hh"
#include "G4ReactionProduct.hh"

class G4NeutronHPNBodyPhaseSpace : public G4VNeutronHPEnergyAngular
{
  public:
  
  G4NeutronHPNBodyPhaseSpace(){}
  ~G4NeutronHPNBodyPhaseSpace(){}
  
  public:
  
  void Init(G4double aMass, G4int aCount)
  {
    theTotalMass=aMass;
    theTotalCount=aCount;
  }

  void Init(ifstream & aDataFile)
  {
    aDataFile >> theTotalMass >> theTotalCount;
    theTotalMass *= G4Neutron::Neutron()->GetPDGMass();
  }
   
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
    
  private:
  
  inline G4double Prob(G4double anEnergy, G4double eMax, G4int n)
  {
    G4double result;
    result = sqrt(anEnergy)*pow(eMax-anEnergy, 3.*n/2.-4.);
    return result;
  }
  
  inline G4double C(G4double anEnergy, G4double mass)
  {
    G4double result;
    if(theTotalCount==3) result = 4./pi/pow(GetEmax(anEnergy, mass),2);
    if(theTotalCount==4) result = 105./32./pow(GetEmax(anEnergy, mass), 3.5);
    if(theTotalCount==5) result = 256./14./pi/pow(GetEmax(anEnergy, mass), 5.);
    return result;
  }
  
  inline G4double GetEmax(G4double anEnergy, G4double mass)
  {
    G4double result;
    G4double tMass = GetTarget()->GetMass();
    G4double pMass = GetNeutron()->GetMass();
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
