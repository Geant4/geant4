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
// $Id: G4NeutronHPProduct.hh,v 1.11 2006/06/29 20:49:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
#ifndef G4NeutronHPProduct_h
#define G4NeutronHPProduct_h 1

#include "G4HadronicException.hh"
#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
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
  
  inline void Init(std::ifstream & aDataFile)
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
      throw G4HadronicException(__FILE__, __LINE__, "distribution law unknown to G4NeutronHPProduct");
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
