// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPProduct.hh,v 1.3 1999-07-02 09:59:54 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPProduct_h
#define G4NeutronHPProduct_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"
#include "G4VNeutronHPEnergyAngular.hh"
#include "G4ReactionProductVector.hh"

#include "G4NeutronHPContEnergyAngular.hh"
#include "G4NeutronHPDiscreteTwoBody.hh"
#include "G4NeutronHPIsotropic.hh"
#include "G4NeutronHPNBodyPhaseSpace.hh"
#include "G4NeutronHPLabAngularEnergy.hh"

class G4NeutronHPProduct
{
  public:
  G4NeutronHPProduct()
  {
    theDist = NULL;
  }
  ~G4NeutronHPProduct()
  {
    if(theDist != NULL) delete theDist;
  }
  
  inline void Init(ifstream & aDataFile)
  {
    aDataFile >> theMassCode>>theMass>>theIsomerFlag>>theDistLaw
              >> theGroundStateQValue>>theActualStateQValue;
    theGroundStateQValue*= eV;
    theActualStateQValue*= eV;
    theYield.Init(aDataFile, eV);
    if(theDistLaw==0)
    {
      // distribution not known, use E-independent, isotropic angular distribution
      theDist = new G4NeutronHPIsotropic;
    }
    else if(theDistLaw == 1)
    {
      // Continuum energy-angular distribution
      theDist = new G4NeutronHPContEnergyAngular;
    }
    else if(theDistLaw == 2)
    {
      // Discrete 2-body scattering
      theDist = new G4NeutronHPDiscreteTwoBody;
    }
    else if(theDistLaw == 3)
    {
      // Isotropic emission
      theDist = new G4NeutronHPIsotropic;
    }
    else if(theDistLaw == 4)
    {
      // Discrete 2-body recoil modification
      // not used for now. @@@@
      theDist = new G4NeutronHPDiscreteTwoBody; 
      // the above is only temporary;
      // recoils need to be addressed
      // properly
      delete theDist;
      theDist = NULL;
    }
    else if(theDistLaw == 5)
    {
      // charged particles only, to be used in a later stage. @@@@
    }
    else if(theDistLaw == 6)
    {
      // N-Body phase space
      theDist = new G4NeutronHPNBodyPhaseSpace;
    }
    else if(theDistLaw == 7)
    {
      // Laboratory angular energy paraetrisation
      theDist = new G4NeutronHPLabAngularEnergy;
    }
    else
    {
      G4Exception("distribution law unknown to G4NeutronHPProduct");
    }
    if(theDist!=NULL)
    {
      theDist->SetQValue(theActualStateQValue);      
      theDist->Init(aDataFile);
    }
  }
  
  G4ReactionProductVector * Sample(G4double anEnergy);
  
  G4double GetMeanYield(G4double anEnergy)
  {
    return theYield.GetY(anEnergy);
  }
  
  void SetNeutron(G4ReactionProduct * aNeutron) 
  { 
    theNeutron = aNeutron; 
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  { 
    theTarget = aTarget; 
  }
  
  inline G4ReactionProduct * GetTarget() { return theTarget; }
  
  inline G4ReactionProduct * GetNeutron() { return theNeutron; }
  
  inline G4double MeanEnergyOfThisInteraction() 
  { 
    G4double result;
    if(theDist == NULL)
    {
      result = 0;
    }
    else
    {
      result=theDist->MeanEnergyOfThisInteraction();
      result *= theCurrentMultiplicity;
    }
    return result;
  }
  
  inline G4double GetQValue() { return theActualStateQValue; }
  private:
   
   // data members

   G4double theMassCode;
   G4double theMass;
   G4int theIsomerFlag;
   G4double theGroundStateQValue;
   G4double theActualStateQValue;
   G4int theDistLaw;  // redundant
   G4NeutronHPVector theYield;
   G4VNeutronHPEnergyAngular *  theDist;
   
   // Utility quantities
   
   G4ReactionProduct * theTarget;
   G4ReactionProduct * theNeutron;

   // cashed values
   
   G4int theCurrentMultiplicity;
  
};

#endif
