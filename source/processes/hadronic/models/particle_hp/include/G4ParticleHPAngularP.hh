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
#ifndef G4ParticleHPAngularP_h
#define G4ParticleHPAngularP_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4InterpolationManager.hh"
#include "G4ParticleHPInterpolator.hh"

class G4ParticleHPAngularP
{
  public:
  
  G4ParticleHPAngularP()
  {
    theCosTh = NULL;
    theProb = NULL;
    theEnergy = 0.;
    nCoeff = 0;
  }
  ~G4ParticleHPAngularP()
  {
    if(theCosTh!=0) delete [] theCosTh;
    if(theProb!=0) delete [] theProb;
  }
  
  inline void Init(std::istream & aDataFile)
  {
    G4double eNeu, cosTheta, probDist;
    G4int nProb;
    aDataFile >> eNeu >> nProb;
    theManager.Init(aDataFile);
    eNeu *= CLHEP::eV;
    Init(eNeu, nProb);
    for (G4int iii=0; iii<nProb; iii++)
    {
      aDataFile >> cosTheta >> probDist;
      SetCosTh(iii, cosTheta);
      SetProb(iii,probDist);
    }  
  }
  
  inline void Init(G4double e, G4int n)
  {
    theCosTh = new G4double[n];
    theProb = new G4double[n];
    theEnergy = e;
    nCoeff = n;
  }
  
  inline void SetEnergy(G4double energy){ theEnergy = energy; }
  inline void SetCosTh(G4int l, G4double coeff) {theCosTh[l]=coeff;}
  inline void SetProb(G4int l, G4double coeff) {theProb[l]=coeff;}
  
  inline G4double GetCosTh(G4int l) {return theCosTh[l];}
  inline G4double GetProb(G4int l) {return theProb[l];}
  inline G4double GetEnergy(){return theEnergy;}
  inline G4int GetNumberOfPoints(){ return nCoeff; }
  inline G4double GetCosTh()
  {
    G4int i;
    G4double rand = G4UniformRand();
    G4double run=0, runo=0;
    for (i=0; i<GetNumberOfPoints(); i++)
    {
      runo=run;
      run += GetProb(i);
      if(run>rand) break;
    }
    if(i == GetNumberOfPoints()) i--;
    G4double costh = theInt.Interpolate(theManager.GetScheme(i), rand, 
                                        runo, run, GetCosTh(i-1), GetCosTh(i));
    return costh;
  }
  
  private:
  
  G4double theEnergy; // neutron energy
  G4ParticleHPInterpolator theInt; // knows tointerpolate
  G4int nCoeff;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4double * theCosTh;
  G4double * theProb; // probability distribution as fcn of theta
                      // integral normalised to 1.
};
#endif
