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
#ifndef G4ParticleHPLegendreTable_h
#define G4ParticleHPLegendreTable_h 1

#include <fstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4InterpolationManager.hh"

class G4ParticleHPLegendreTable
{
  public:
  G4ParticleHPLegendreTable()
  {
    nCoeff=0; 
    theCoeff = 0;
    theRep = 0;
    theEnergy = 0.0;
    theTemp = 0.0;
  }
  ~G4ParticleHPLegendreTable(){if(theCoeff!=0) delete [] theCoeff;}
  
  void operator= (const G4ParticleHPLegendreTable & aSet)
  {
    if(&aSet!=this)
    {
      theRep = aSet.theRep;
      theEnergy = aSet.theEnergy;
      theTemp = aSet.theTemp;
      theManager = aSet.theManager;
      nCoeff = aSet.nCoeff;
      if(theCoeff!=0) delete [] theCoeff;
      theCoeff = new G4double[nCoeff];
      for(G4int i=0; i<nCoeff; i++)
      {
        theCoeff[i] = aSet.theCoeff[i];
      }
    }
  }
  
  inline void Init(std::istream & aDataFile) 
  {
    G4double eNeu, coeff;
    G4int nPoly;
    aDataFile >> eNeu >> nPoly;
    eNeu *= CLHEP::eV;
    Init(eNeu, nPoly);
    for(G4int l=0; l<nPoly; l++)
    {
      aDataFile >> coeff;
      SetCoeff(l+1, coeff);
    }
  }
  
  inline void Init(G4double e, G4int n)
  {
    nCoeff = n+1;
    theCoeff = new G4double[nCoeff];
    for(G4int i=0; i<nCoeff; i++) theCoeff[i] = 0;
    theCoeff[0]=1.;
    theEnergy = e;
//    G4cout << "G4ParticleHPLegendreTable::Init called "<<e<<" "<<n<<G4endl;
  }
  inline void SetEnergy(G4double energy){ theEnergy = energy; }
  inline void SetTemperature(G4double temp){ theTemp = temp; }
  inline void SetCoeff(G4int l, G4double coeff) {theCoeff[l]=coeff;}
  inline void SetRepresentation(G4int aRep) {theRep = aRep;}
  
  inline G4double GetCoeff(G4int l) {return theCoeff[l];}
  inline G4double GetEnergy(){return theEnergy;}
  inline G4double GetTemperature(){return theTemp;}
  inline G4int GetNumberOfPoly() {return nCoeff;}
  inline G4int GetRepresentation() {return theRep;}
  inline const G4InterpolationManager & GetManager() {return theManager;}
  private:
  
  G4int theRep;
  G4double theEnergy;
  G4double theTemp;
  G4int nCoeff;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4double * theCoeff;
};

#endif
