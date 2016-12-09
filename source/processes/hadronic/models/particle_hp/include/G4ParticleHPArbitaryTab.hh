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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPArbitaryTab_h
#define G4ParticleHPArbitaryTab_h 1

#include <fstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4ParticleHPVector.hh"
#include "G4VParticleHPEDis.hh"
#include "G4InterpolationManager.hh"

// we will need a List of these .... one per term.

class G4ParticleHPArbitaryTab : public G4VParticleHPEDis
{
  public:
  G4ParticleHPArbitaryTab()
  {
   theDistFunc = 0;
   nDistFunc = 0;
  }
  ~G4ParticleHPArbitaryTab()
  {
   if(theDistFunc!=0) delete [] theDistFunc;
  }
  
  inline void Init(std::istream & theData)
  {
    G4int i;
    theFractionalProb.Init(theData, CLHEP::eV);
    theData >> nDistFunc; // = number of incoming n energy points
    theDistFunc = new G4ParticleHPVector [nDistFunc];
    theManager.Init(theData);
    G4double currentEnergy;
    for(i=0; i<nDistFunc; i++)
    {
      theData >> currentEnergy;
      theDistFunc[i].SetLabel(currentEnergy*CLHEP::eV);
      theDistFunc[i].Init(theData, CLHEP::eV);
      theDistFunc[i].IntegrateAndNormalise();
      //************************************************************************
      //EMendoza:
      //ThinOut() assumes that the data is linear-linear, what is false:
      //theDistFunc[i].ThinOut(0.02); // @@@ optimization to be finished.
      //************************************************************************
    }

    //************************************************************************
    //EMendoza:
    //Here we calculate the thresholds for the 2D sampling:
    for(i=0; i<nDistFunc; i++){
      G4int np=theDistFunc[i].GetVectorLength();
      theLowThreshold[i]=theDistFunc[i].GetEnergy(0);
      theHighThreshold[i]=theDistFunc[i].GetEnergy(np-1);
      for(G4int j=0;j<np-1;j++){
	if(theDistFunc[i].GetXsec(j+1)>1.e-20){
	  theLowThreshold[i]=theDistFunc[i].GetEnergy(j);
	  break;
	}
      }
      for(G4int j=1;j<np;j++){
	if(theDistFunc[i].GetXsec(j-1)>1.e-20){
	  theHighThreshold[i]=theDistFunc[i].GetEnergy(j);
	}
      }
    }
     //************************************************************************
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy) ;
  
  private:
  
  G4ParticleHPVector theFractionalProb;
  G4int nDistFunc;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4ParticleHPVector * theDistFunc; // one per incoming energy
  G4ParticleHPVector theBuffer;
  //************************************************************************
  //EMendoza:
  G4double theLowThreshold[1000];
  G4double theHighThreshold[1000];
  //************************************************************************

};

#endif
