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
#ifndef ANAParticle_h
#define ANAParticle_h

#include "globals.hh"
#include <fstream>

class ANAParticle
{
  public:
    ANAParticle(G4int aCode, G4double apx, G4double apy, G4double apz, G4double anenergy):
      charge(0), pdgCode(aCode), energy(anenergy), px(apx), py(apy), pz(apz) {}

    ANAParticle(){}

    G4bool Init(ifstream& theData);
    
    G4int GetCharge()           {return charge;}
    G4int GetPDGCode()          {return pdgCode;}
    G4double GetEnergy()        {return energy;}
    G4double GetKineticEnergy() {return energy - GetMomentum();}
    G4double GetPx()            {return px;}
    G4double GetPy()            {return py;}
    G4double GetPz()            {return pz;}
    G4double GetMomentum()      {return sqrt(max(0., px*px+py*py+pz*pz));}
    G4double GetCosTheta() {G4double p=GetMomentum(); if(p)return pz/p; return 0.;}
    G4double GetWeight() {G4double p=GetMomentum(); if(p)return 1./(twopi*p); return 0.;}
    
  private:
  G4int charge;
  G4int pdgCode;
  G4double energy;
  G4double px;
  G4double py;
  G4double pz;
};

inline G4bool ANAParticle::Init(ifstream& theData)
{
  G4int dummy;
  if(!(theData>>charge)) return false;
  theData>>pdgCode;
  theData>>energy;
  theData>>px;
  theData>>py;
  theData>>pz;
  theData>>dummy;
  return true;
}

#endif
