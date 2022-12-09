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
//
// 080718 Add ClearHistories method and related class member
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPContAngularPar_h
#define G4ParticleHPContAngularPar_h 1

#include <fstream>
#include <set>

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleHPList.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4InterpolationManager.hh"
#include "G4Cache.hh"

class G4ParticleDefinition;

class G4ParticleHPContAngularPar
{

   struct toBeCached
   {
      G4bool fresh;
      G4double currentMeanEnergy;
      G4double remaining_energy; 
      G4double theTargetCode;
      G4ReactionProduct* theTarget;
      G4ReactionProduct* thePrimary;
      toBeCached()
      : fresh(true),currentMeanEnergy(-2.0),remaining_energy(0.0),
        theTargetCode(-1.0),theTarget(0),thePrimary(0) {}
   };

  public:
  
  G4ParticleHPContAngularPar()
  {
    theAngular = 0;
    //currentMeanEnergy = -2;
    //fresh = true;
    fCache.Put(0);
    theMinEner = DBL_MAX;
    theMaxEner = -DBL_MAX;
    theEnergy = -1;
    nEnergies = -1;
    nDiscreteEnergies = -1;
    nAngularParameters = -1;
    theProjectile = 0;
    adjustResult = true;
  }

  G4ParticleHPContAngularPar(G4ParticleHPContAngularPar & val)
  {
    theEnergy         = val.theEnergy;
    nEnergies         = val.nEnergies;
    nDiscreteEnergies = val.nDiscreteEnergies;
    nAngularParameters= val.nAngularParameters;
    theProjectile     = val.theProjectile;
    theManager        = val.theManager;
    theInt            = val.theInt;
    adjustResult      = val.adjustResult;
    theMinEner        = val.theMinEner;
    theMaxEner        = val.theMaxEner;
    theEnergiesTransformed = val.theEnergiesTransformed;
    theDiscreteEnergies = val.theDiscreteEnergies;
    theDiscreteEnergiesOwn = val.theDiscreteEnergiesOwn;
    fCache.Put(0);
    theAngular        = new G4ParticleHPList[nEnergies];
    for(G4int ie=0;ie<nEnergies;++ie) {
      theAngular[ie].SetLabel(val.theAngular[ie].GetLabel());
      for(G4int ip=0;ip<nAngularParameters;++ip) {
	theAngular[ie].SetValue(ip,val.theAngular[ie].GetValue(ip));
      }
    }
  }

  G4ParticleHPContAngularPar(G4ParticleDefinition* projectile); 

  ~G4ParticleHPContAngularPar()
  {
    if (theAngular !=0 ) delete [] theAngular;
    if (fCache.Get() != 0) delete fCache.Get();
  }
  
  void Init(std::istream & aDataFile, G4ParticleDefinition* projectile);
  
  G4ReactionProduct* Sample(G4double anEnergy, G4double massCode, G4double mass, 
                            G4int angularRep, G4int interpol);
  
  G4double GetEnergy() const
  { 
    return theEnergy;
  }
  
  void SetPrimary(G4ReactionProduct * aPrimary)
  {
    fCache.Get()->thePrimary = aPrimary;
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  {
    fCache.Get()->theTarget = aTarget;
  }
  
  void SetTargetCode(G4double aTargetCode)
  {
    fCache.Get()->theTargetCode = aTargetCode;
  }
  
  void SetInterpolation(G4int theInterpolation)
  {
    theManager.Init(theInterpolation, nEnergies); // one range only
  }

  void BuildByInterpolation(G4double anEnergy, G4InterpolationScheme aScheme, 
             G4ParticleHPContAngularPar & store1, 
             G4ParticleHPContAngularPar & store2);
    // NOTE: this interpolates legendre coefficients

  void PrepareTableInterpolation();
  
  G4double MeanEnergyOfThisInteraction()
  {
    G4double result;
    if(fCache.Get()->currentMeanEnergy<-1)
    {
      return 0;
      // throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar: Logical error in Product class");
    }
    else
    {
      result = fCache.Get()->currentMeanEnergy;
    }
    fCache.Get()->currentMeanEnergy = -2;
    return result;
  }
  
  G4int GetNEnergies() const 
  {
    return nEnergies; 
  }
  G4int GetNDiscreteEnergies() const 
  {
    return nDiscreteEnergies; 
  }
  std::set<G4double> GetEnergiesTransformed() const 
  {
    return theEnergiesTransformed;
  }
  G4int GetNEnergiesTransformed() const 
  {
    return (G4int)theEnergiesTransformed.size();
  }
  G4double GetMinEner() const 
  {
    return theMinEner;
  }
  G4double GetMaxEner() const 
  { 
    return theMaxEner;
  }
  std::map<G4double,G4int> GetDiscreteEnergiesOwn() const 
  {
    return theDiscreteEnergiesOwn;
  }
  G4ParticleHPList* GetAngDataList() const
  {
    return theAngular; 
  }
  
  void ClearHistories()
  { 
    if ( fCache.Get() == 0 ) cacheInit();
    fCache.Get()->fresh = true;
  }

  void Dump() const;

private:
  
  // incoming particle
  G4double theEnergy;
  // number of exit channel energies
  G4int nEnergies; 
  // number of discrete exit channels
  G4int nDiscreteEnergies;
  // number of angular paramerers per channel
  G4int nAngularParameters;
  // knows the interpolation between List labels
  G4InterpolationManager theManager; 
  // on per exit-channel energy
  G4ParticleHPList * theAngular; 
  
  G4ParticleHPInterpolator theInt;
  
private:

  G4Cache< toBeCached* > fCache;
  void cacheInit()
  {
    toBeCached* val = new toBeCached;
    val->currentMeanEnergy = -2;
    val->remaining_energy = 0;
    val->fresh=true;
    fCache.Put( val );
  };

  G4ParticleDefinition* theProjectile;

  G4bool adjustResult;
    // if not set it will not force the conservation of energy in angularRep==1,
    // but will sample the particle energy according to the database

  G4double theMinEner;
  G4double theMaxEner;
  std::set<G4double> theEnergiesTransformed;
  std::set<G4double> theDiscreteEnergies;
  std::map<G4double,G4int> theDiscreteEnergiesOwn;
};

#endif
