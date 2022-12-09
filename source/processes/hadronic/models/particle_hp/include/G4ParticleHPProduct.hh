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
#ifndef G4ParticleHPProduct_h
#define G4ParticleHPProduct_h 1

#include "G4HadronicException.hh"
#include "globals.hh"
#include "G4ParticleHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ReactionProductVector.hh"

#include "G4ParticleHPContEnergyAngular.hh"
#include "G4ParticleHPDiscreteTwoBody.hh"
#include "G4ParticleHPIsotropic.hh"
#include "G4ParticleHPNBodyPhaseSpace.hh"
#include "G4ParticleHPLabAngularEnergy.hh"
#include "G4Cache.hh"
class G4ParticleDefinition;

enum G4HPMultiMethod { G4HPMultiPoisson, G4HPMultiBetweenInts };

class G4ParticleHPProduct
{
   struct toBeCached
   {
      G4ReactionProduct* theProjectileRP;
      G4ReactionProduct* theTarget;
      G4int theCurrentMultiplicity;
      toBeCached()
        : theProjectileRP(0), theTarget(0), theCurrentMultiplicity(-1) {}
   };

public:

  G4ParticleHPProduct()
  {
    theDist = 0;
    toBeCached val;
    fCache.Put( val );

    char * method = std::getenv( "G4PHP_MULTIPLICITY_METHOD" );
    if( method )
    {
      if( G4String(method) == "Poisson" ) {
	theMultiplicityMethod = G4HPMultiPoisson;
      } else if( G4String(method) == "BetweenInts" ) {
	theMultiplicityMethod = G4HPMultiBetweenInts;
      } else {
	throw G4HadronicException(__FILE__, __LINE__, ("multiplicity method unknown to G4ParticleHPProduct" + G4String(method)).c_str());
      }
    }
    else
    {
      theMultiplicityMethod = G4HPMultiPoisson;
    }
    theMassCode = 0.0;
    theMass = 0.0;
    theIsomerFlag = 0;
    theGroundStateQValue = 0.0;
    theActualStateQValue = 0.0;
    theDistLaw = -1;
  }

  ~G4ParticleHPProduct()
  { 
    if(theDist != 0) delete theDist;
  }
  
  inline void Init(std::istream & aDataFile, G4ParticleDefinition* projectile)
  {
    aDataFile >> theMassCode>>theMass>>theIsomerFlag>>theDistLaw
              >> theGroundStateQValue>>theActualStateQValue;
    theGroundStateQValue*= CLHEP::eV;
    theActualStateQValue*= CLHEP::eV;
    theYield.Init(aDataFile, CLHEP::eV);
    theYield.Hash();
    if(theDistLaw==0)
    {
      // distribution not known, use E-independent, isotropic
      // angular distribution
      theDist = new G4ParticleHPIsotropic;
    }
    else if(theDistLaw == 1)
    {
      // Continuum energy-angular distribution
      theDist = new G4ParticleHPContEnergyAngular(projectile);
    }
    else if(theDistLaw == 2)
    {
      // Discrete 2-body scattering
      theDist = new G4ParticleHPDiscreteTwoBody;
    }
    else if(theDistLaw == 3)
    {
      // Isotropic emission
      theDist = new G4ParticleHPIsotropic;
    }
    else if(theDistLaw == 4)
    {
      // Discrete 2-body recoil modification
      // not used for now. @@@@
      theDist = new G4ParticleHPDiscreteTwoBody; 
      // the above is only temporary;
      // recoils need to be addressed
      // properly
      delete theDist;
      theDist = 0;
    }
    //    else if(theDistLaw == 5)
    //    {
      // charged particles only, to be used in a later stage. @@@@
    //    }
    else if(theDistLaw == 6)
    {
      // N-Body phase space
      theDist = new G4ParticleHPNBodyPhaseSpace;
    }
    else if(theDistLaw == 7)
    {
      // Laboratory angular energy paraetrisation
      theDist = new G4ParticleHPLabAngularEnergy;
    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "distribution law unknown to G4ParticleHPProduct");
    }
    if(theDist!=0)
    {
      theDist->SetQValue(theActualStateQValue);      
      theDist->Init(aDataFile);
    }
  }
  
  G4int GetMultiplicity(G4double anEnergy);
  G4ReactionProductVector * Sample(G4double anEnergy, G4int nParticles);
  
  G4double GetMeanYield(G4double anEnergy)
  {
    return theYield.GetY(anEnergy);
  }
  
  void SetProjectileRP(G4ReactionProduct * aIncidentPart) 
  { 
    fCache.Get().theProjectileRP = aIncidentPart; 
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  { 
    fCache.Get().theTarget = aTarget; 
  }
  
  inline G4ReactionProduct * GetTarget()
  {
    return fCache.Get().theTarget;
  }
  
  inline G4ReactionProduct * GetProjectileRP()
  {
    return fCache.Get().theProjectileRP;
  }
  
  inline G4double MeanEnergyOfThisInteraction() 
  { 
    G4double result;
    if(theDist == 0)
    {
      result = 0;
    }
    else
    {
      result=theDist->MeanEnergyOfThisInteraction();
      result *= fCache.Get().theCurrentMultiplicity;
    }
    return result;
  }
  
  inline G4double GetQValue()
  {
    return theActualStateQValue;
  }

  //TK120515 For migration of frameFlag (MF6 LCT) = 3 in
  //G4ParticleHPEnAngCorrelation
  G4double GetMassCode() {return theMassCode;}
  G4double GetMass() {return theMass;}

private:
   
  // data members

  G4double theMassCode;
  G4double theMass;
  G4int theIsomerFlag;
  G4double theGroundStateQValue;
  G4double theActualStateQValue;
  G4int theDistLaw;  // redundant
  G4ParticleHPVector theYield;
  G4VParticleHPEnergyAngular *  theDist;
   
  // cashed values
  //
  G4Cache<toBeCached> fCache;

  G4HPMultiMethod theMultiplicityMethod;  
};

#endif
