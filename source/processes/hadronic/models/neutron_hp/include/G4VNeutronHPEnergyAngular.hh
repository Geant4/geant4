// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VNeutronHPEnergyAngular.hh,v 1.4 1999-12-15 14:53:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4VNeutronHPEnergyAngular_h
#define G4VNeutronHPEnergyAngular_h 1

#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4ReactionProduct.hh"

class G4VNeutronHPEnergyAngular
{
  public:
  
  G4VNeutronHPEnergyAngular()
  {
    theTarget = NULL;
    theNeutron = NULL;
    theQValue=0;
  }
  virtual ~G4VNeutronHPEnergyAngular(){}
  
  public:
  
  virtual void Init(G4std::ifstream & aDataFile) = 0;
  virtual G4ReactionProduct * Sample(G4double anEnergy, 
                                     G4double massCode, 
                                     G4double mass) = 0;
  virtual G4double MeanEnergyOfThisInteraction() = 0; // returns value cashed in sample
  
  void SetNeutron(G4ReactionProduct * aNeutron) 
  { 
    theNeutron = aNeutron; 
    if(theTarget!=NULL) theCMS = *theNeutron+*theTarget;
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  { 
    theTarget = aTarget; 
  }
  
  G4ReactionProduct * GetTarget() { return theTarget; }
  
  G4ReactionProduct * GetNeutron() { return theNeutron; }
  
  G4ReactionProduct * GetCMS() { return &theCMS; }

  inline void SetQValue(G4double aValue) { theQValue = aValue; }
  
  protected:
  
  inline G4double GetQValue() { return theQValue; }
  
  private:
  
  G4double theQValue;
    
  G4ReactionProduct * theTarget;
  G4ReactionProduct * theNeutron;
  G4ReactionProduct theCMS;
    
};
#endif
