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
// 080721 Add ClearHistories() method by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPContEnergyAngular_h
#define G4ParticleHPContEnergyAngular_h 1

#include <fstream>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ParticleHPContAngularPar.hh"
#include "G4InterpolationManager.hh"
#include "G4Cache.hh"

class G4ParticleDefinition;

// we will need one of these per product.

class G4ParticleHPContEnergyAngular : public G4VParticleHPEnergyAngular
{
public:
  
  G4ParticleHPContEnergyAngular(G4ParticleDefinition* proj)
    : theProjectile(proj)
  {
    theAngular = 0;
    currentMeanEnergy.Put( -2 );
    theTargetCode = -1.0;
    theAngularRep = -1;
    nEnergy = -1;
    theInterpolation = -1;
    fCacheAngular.Put(0); //fix
  }
  
  ~G4ParticleHPContEnergyAngular()
  {
    if(theAngular!=0) delete [] theAngular;
    if (fCacheAngular.Get() != 0) delete fCacheAngular.Get(); //fix
  }
  
  void Init(std::istream & aDataFile)
  {
    aDataFile >> theTargetCode >> theAngularRep >> theInterpolation >> nEnergy;
    theAngular = new G4ParticleHPContAngularPar[nEnergy];
    theManager.Init(aDataFile);
    for(G4int i=0; i<nEnergy; i++)
    {
      theAngular[i].Init(aDataFile, theProjectile);
      theAngular[i].SetInterpolation(theInterpolation);
#ifndef PHP_AS_HP
      theAngular[i].PrepareTableInterpolation();
#endif
    }
  }

  G4double MeanEnergyOfThisInteraction();
  G4ReactionProduct* Sample(G4double anEnergy, G4double massCode, G4double mass);
  void ClearHistories(); 
  
private:
  
  G4double theTargetCode;
  G4int theAngularRep;
  G4int nEnergy;
  
  G4int theInterpolation;

  G4InterpolationManager theManager; // knows the interpolation between stores
  G4ParticleHPContAngularPar * theAngular;
  
  G4Cache<G4double> currentMeanEnergy;
  G4Cache<G4ParticleHPContAngularPar*> fCacheAngular; //fix
  G4ParticleDefinition* theProjectile;
};
#endif
