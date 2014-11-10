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
// 080718 Add ClearHistories() method by T. Koi
//
#ifndef G4VNeutronHPEnergyAngular_h
#define G4VNeutronHPEnergyAngular_h 1

#include <fstream>

#include "globals.hh"
#include "G4ios.hh"
#include "G4ReactionProduct.hh"

#include "G4Cache.hh"

class G4VNeutronHPEnergyAngular
{

   struct toBeCached {
      G4ReactionProduct* theNeutron;
      G4ReactionProduct* theTarget;
      G4ReactionProduct* theCMS;
      toBeCached() : theNeutron(NULL),theTarget(NULL),theCMS(NULL) {};
   };

  public:
  
  G4VNeutronHPEnergyAngular()
  {
    //theTarget = 0;
    //theNeutron = 0;
    theQValue=0;
      toBeCached val;
      fCache.Put( val );
  }
  virtual ~G4VNeutronHPEnergyAngular(){}
  
  public:
  
  virtual void Init(std::istream & aDataFile) = 0;
  virtual G4ReactionProduct * Sample(G4double anEnergy, 
                                     G4double massCode, 
                                     G4double mass) = 0;
  virtual G4double MeanEnergyOfThisInteraction() = 0; // returns value cashed in sample
  
  void SetNeutron(G4ReactionProduct * aNeutron) 
  { 
    fCache.Get().theNeutron = aNeutron; 
    //if(theTarget!=0) theCMS = *theNeutron+*theTarget;
    //if ( fCache.Get().theTarget != 0 ) theCMS = *fCache.Get().theNeutron+*fCache.Get().theTarget;
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  { 
    fCache.Get().theTarget = aTarget; 
  }
  
  G4ReactionProduct * GetTarget() { return fCache.Get().theTarget; }
  
  G4ReactionProduct * GetNeutron() { return fCache.Get().theNeutron; }
  
  //G4ReactionProduct * GetCMS() { return &theCMS; }
  G4ReactionProduct * GetCMS() { 
     *fCache.Get().theCMS = *fCache.Get().theNeutron + *fCache.Get().theTarget;
     return fCache.Get().theCMS; 
  }

  inline void SetQValue(G4double aValue) { theQValue = aValue; }
  
  protected:
  
  inline G4double GetQValue() { return theQValue; }
  
  private:
  
  G4double theQValue;
    
  //G4ReactionProduct * theTarget;
  //G4ReactionProduct * theNeutron;
  //G4ReactionProduct theCMS;

   private:
      G4Cache<toBeCached> fCache;
    
   public:
      virtual void ClearHistories(){;};
};
#endif
