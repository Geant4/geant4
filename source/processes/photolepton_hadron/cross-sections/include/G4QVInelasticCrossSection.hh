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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QVInelasticCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//

#ifndef G4QVInelasticCrossSection_h
#define G4QVInelasticCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>

class G4QVInelasticCrossSection
{
public:

  G4QVInelasticCrossSection()  {}
  virtual ~G4QVInelasticCrossSection() {}

  virtual G4double GetCrossSection(G4double Energy, G4int Z, G4int N) = 0;
  virtual G4double ThresholdEnergy(G4int Z, G4int N) = 0; // #0 only for PhotoNuclearCS
  virtual G4int GetFunctions(G4double A, G4double* y, G4double* z) = 0; // Tablated CS
  // Linear fit for YN[N] tabulated (from X0 with fixed step DX) function to X point
  G4double EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, const G4double* Y)
  {
    if(DX<=0.||N<2){G4cout<<"*G4QVInelCS::EquLinFit:D="<<DX<<",N="<<N<<G4endl;return Y[0];}
    G4double d=(X-X0)/DX; G4int j=static_cast<int>(d); // Calculate numberOfChannels in Tab
    G4int  N2=N-2; if(j<0) j=0; else if(j>N2) j=N2;    // Restrict theSelectedChannelNumber
    d-=j; G4double yi=Y[j]; G4double s=yi+(Y[j+1]-yi)*d; // Calculate value of the sigma
    if(s<0.)G4cout<<"G4QInCS::EquLinFit:j="<<j<<",i="<<yi<<",a="<<Y[j+1]<<",d="<<d<<G4endl;
    return s;
  }

};

#endif
