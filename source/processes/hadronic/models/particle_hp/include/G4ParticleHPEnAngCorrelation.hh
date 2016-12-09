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
#ifndef G4ParticleHPEnAngCorrelation_h
#define G4ParticleHPEnAngCorrelation_h 1

#include "globals.hh"
#include "G4ParticleHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4ParticleHPProduct.hh"
#include "G4ReactionProduct.hh"
#include "G4Cache.hh"
class G4ParticleDefinition;

class G4ParticleHPEnAngCorrelation
{

   struct toBeCached {
      G4ReactionProduct* theProjectileRP;
      G4ReactionProduct* theTarget;
      G4double theTotalMeanEnergy;
      toBeCached() : theProjectileRP(NULL),theTarget(NULL),theTotalMeanEnergy(-1.0) {};
   };

  public:
  G4ParticleHPEnAngCorrelation() // for G4ParticleHPCaptureFS::theMF6FinalState
  {
    theProjectile = G4Neutron::Neutron();
    theProducts = 0;
    inCharge = false;
    toBeCached val;
    fCache.Put( val );
    //theTotalMeanEnergy = -1.;
    fCache.Get().theTotalMeanEnergy = -1.;
    targetMass = 0.0;
    frameFlag = 0;
    nProducts = 0;
    bAdjustFinalState = true;
  }
  G4ParticleHPEnAngCorrelation(G4ParticleDefinition* proj)
    : theProjectile(proj)
  {
    theProducts = 0;
    inCharge = false;
    toBeCached val;
    fCache.Put( val );
    //theTotalMeanEnergy = -1.;
    fCache.Get().theTotalMeanEnergy = -1.;
    targetMass = 0.0;
    frameFlag = 0;
    nProducts = 0;
    bAdjustFinalState = true;
  }

  ~G4ParticleHPEnAngCorrelation()
  {
    if(theProducts!=0) delete [] theProducts;
  }
  
  inline void Init(std::istream & aDataFile)
  {
    bAdjustFinalState = true;
    const char* ctmp = getenv("G4PHP_DO_NOT_ADJUST_FINAL_STATE");
    if( ctmp && G4String(ctmp) == "1" ) {
      bAdjustFinalState = false;
    }
//T.K. Comment out following line to keep the condition at the  validation efforts compairng NeutronHP and PartileHP for neutrons (2015 Sep.)
//#ifdef PHP_AS_HP
//    bAdjustFinalState = false;
//#endif
 
    inCharge = true;
    aDataFile>>targetMass>>frameFlag>>nProducts;
    theProducts = new G4ParticleHPProduct[nProducts];
    for(G4int i=0; i<nProducts; i++)
    {
      theProducts[i].Init(aDataFile,theProjectile);
    }

  }
  
  G4ReactionProduct * SampleOne(G4double anEnergy);

  G4ReactionProductVector * Sample(G4double anEnergy);
  
  inline void SetTarget(G4ReactionProduct & aTarget)
  {
    fCache.Get().theTarget = &aTarget;
    for(G4int i=0;i<nProducts;i++)theProducts[i].SetTarget(fCache.Get().theTarget);
  }
  
  inline void SetProjectileRP(G4ReactionProduct & aIncidentPart)
  {
    fCache.Get().theProjectileRP = &aIncidentPart;
    for(G4int i=0;i<nProducts;i++)theProducts[i].SetProjectileRP(fCache.Get().theProjectileRP);
  }
  
  inline G4bool InCharge()
  {
    return inCharge;
  }
  
  inline G4double GetTargetMass() { return targetMass; }
  
  G4double GetTotalMeanEnergy()
  {
     // cashed in 'sample' call
    return fCache.Get().theTotalMeanEnergy; 
  }
private:
   
  // data members
  
  G4double targetMass;
  G4int frameFlag; // =1: Target rest frame; =2: CMS system; incident always in lab
  G4int nProducts;
  G4ParticleHPProduct * theProducts;
  G4bool inCharge;
    
  // Utility quantities
  
  //G4ReactionProduct theTarget;
  //G4ReactionProduct theProjectileRP;
  
  // cashed values
  
  //G4double theTotalMeanEnergy;
     G4Cache<toBeCached> fCache;

  G4ParticleDefinition* theProjectile;

  G4bool bAdjustFinalState;  

};

#endif
