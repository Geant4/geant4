//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPAngularP.hh,v 1.7 2002-12-12 19:18:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPAngularP_h
#define G4NeutronHPAngularP_h 1

#include "globals.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPAngularP
{
  public:
  
  G4NeutronHPAngularP()
  {
    theCosTh = NULL;
    theProb = NULL;
  }
  ~G4NeutronHPAngularP()
  {
    if(theCosTh!=NULL) delete [] theCosTh;
    if(theProb!=NULL) delete [] theProb;
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    G4double eNeu, cosTheta, probDist;
    G4int nProb;
    aDataFile >> eNeu >> nProb;
    theManager.Init(aDataFile);
    eNeu *= eV;
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
  G4NeutronHPInterpolator theInt; // knows tointerpolate
  G4int nCoeff;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4double * theCosTh;
  G4double * theProb; // probability distribution as fcn of theta
                      // integral normalised to 1.
};
#endif
