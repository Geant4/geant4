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
// 070613 fix memory leaking by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPAngular_h
#define G4ParticleHPAngular_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4ReactionProduct.hh"
#include "Randomize.hh"
#include "G4ParticleHPLegendreStore.hh"
#include "G4ParticleHPPartial.hh"
#include "G4Cache.hh"

class G4ParticleHPAngular
{

   struct toBeCached {
      const G4ReactionProduct* theProjectileRP;
      const G4ReactionProduct* theTarget;
      toBeCached() : theProjectileRP(NULL),theTarget(NULL) {};
   };

    public:
    
  G4ParticleHPAngular()
  {
//TKDB110511
    //theAngularDistributionType = 0;
    //theIsoFlag = false;
    theIsoFlag = true;
// TKDB
      theCoefficients = 0;
      theProbArray = 0;
      toBeCached val;
      fCache.Put( val );
      
      theAngularDistributionType = 0;
      frameFlag = 0;
      targetMass = 0.0;
  } 

  ~G4ParticleHPAngular()
   {
// TKDB
      delete theCoefficients;
      delete theProbArray;
   }
  
  void Init(std::istream & aDataFile);
  
  void SampleAndUpdate(G4ReactionProduct & anIncidentParticle);
    
  void SetTarget(const G4ReactionProduct & aTarget) { fCache.Get().theTarget = &aTarget; }

  void SetProjectileRP(const G4ReactionProduct & anIncidentParticleRP) { fCache.Get().theProjectileRP = &anIncidentParticleRP; }

  inline G4double GetTargetMass() { return targetMass; }

  private:
  
  // the type of distribution; currently 
  // isotropic (0), 
  // and legendre representation (1)
  // probability distribution (2)
  // are supported
  
  G4int theAngularDistributionType;
  G4int frameFlag; // 1=Lab, 2=CMS
    
  G4bool theIsoFlag; // isotropic or not?
  
  G4ParticleHPLegendreStore * theCoefficients; // the legendre coefficients

  G4ParticleHPPartial * theProbArray; // the probability array p,costh for energy

  private:
  
  G4double targetMass;

  //G4ReactionProduct theTarget;
  //G4ReactionProduct theProjectileRP;
     G4Cache<toBeCached> fCache;

};

#endif
