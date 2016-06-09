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
// $Id: G4EmCorrections.hh,v 1.3 2005/02/26 22:01:20 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmCorrections
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.01.2005
//
// Modifications:
//
//
// Class Description:
//
// This class provides calculation of EM corrections to ionisation 
//

// -------------------------------------------------------------------
//

#ifndef G4EmCorrections_h
#define G4EmCorrections_h 1

#include "globals.hh"
#include "G4AtomicShells.hh"
#include "G4ionEffectiveCharge.hh"

class G4Material;
class G4ParticleDefinition;

class G4EmCorrections 
{

public:

  G4EmCorrections();

  virtual ~G4EmCorrections();

  G4double HighOrderCorrections(const G4ParticleDefinition* p,
                                const G4Material* material,
                                      G4double kineticEnergy);

  G4double Bethe(const G4ParticleDefinition* p,
                 const G4Material* material,
                       G4double kineticEnergy);

  G4double SpinCorrection(const G4ParticleDefinition* p,
                          const G4Material* material,
                                G4double kineticEnergy);

  G4double KShellCorrection(const G4ParticleDefinition* p,
                            const G4Material* material,
                                  G4double kineticEnergy);

  G4double LShellCorrection(const G4ParticleDefinition* p,
                            const G4Material* material,
                                  G4double kineticEnergy);

  G4double ShellCorrection(const G4ParticleDefinition* p,
                           const G4Material* material,
                                 G4double kineticEnergy);

  G4double ShellCorrectionSTD(const G4ParticleDefinition* p,
                              const G4Material* material,
                                    G4double kineticEnergy);

  G4double DensityCorrection(const G4ParticleDefinition* p,
                             const G4Material* material,
                                   G4double kineticEnergy);

  G4double BarkasCorrection(const G4ParticleDefinition* p,
                            const G4Material* material,
                                  G4double kineticEnergy);

  G4double BlochCorrection(const G4ParticleDefinition* p,
                           const G4Material* material,
                                 G4double kineticEnergy);

  G4double MottCorrection(const G4ParticleDefinition* p,
                          const G4Material* material,
                                G4double kineticEnergy);

  G4double NuclearDEDX(const G4ParticleDefinition* p,
                       const G4Material* material,
 		             G4double kineticEnergy,
                             G4bool fluct = true);

private:

  void Initialise();

  G4double KShell(G4double theta, G4double eta); 

  G4double LShell(G4double theta, G4double eta); 

  G4int Index(G4double x, G4double* y, G4int n); 

  G4double Value(G4double xv, G4double x1, G4double x2, G4double y1, G4double y2); 

  G4double Value2(G4double xv, G4double yv, G4double x1, G4double x2,
                  G4double y1, G4double y2,
                  G4double z11, G4double z21, G4double z12, G4double z22); 

  G4double NuclearStoppingPower(G4double e, G4double z1, G4double z2, 
                                            G4double m1, G4double m2);

  // hide assignment operator
  G4EmCorrections & operator=(const G4EmCorrections &right);
  G4EmCorrections(const G4EmCorrections&);

  G4double     engBarkas[47];
  G4double     corBarkas[47];
  G4double     e[104];
  G4double     a[104];
  G4double     theZieglerFactor;
  G4double     alpha2;
  G4bool       lossFlucFlag;

  G4int        nK;
  G4int        nL;
  G4int        nEtaK;
  G4int        nEtaL;

  G4double     COSEB[14];
  G4double     COSXI[14];
  G4double     ZD[11];

  G4double     TheK[20];
  G4double     SK[20];
  G4double     TK[20];
  G4double     UK[20];
  G4double     VK[20];
  G4double     ZK[20];

  G4double     TheL[26];
  G4double     SL[26];
  G4double     TL[26];
  G4double     UL[26];
  G4double     VL[26];

  G4double     Eta[29];
  G4double     CK[20][29];
  G4double     CL[26][28];
  G4double     HM[53];
  G4double     HN[31];
  G4double     MSH[93];
  G4double     TAU[93];

  G4AtomicShells        shells;
  G4ionEffectiveCharge  effCharge; 
};

inline G4int G4EmCorrections::Index(G4double x, G4double* y, G4int n) 
{
  G4int idx = n-1;
  do {idx--;} while (idx>0 && x<y[idx]);
  return idx;
}

inline G4double G4EmCorrections::Value(G4double xv, G4double x1, G4double x2, 
                                       G4double y1, G4double y2) 
{
  return y1 + (y2 - y1)*(xv - x1)/(x2 - x1);
}

inline G4double G4EmCorrections::Value2(G4double xv, G4double yv, G4double x1, G4double x2,
                                        G4double y1, G4double y2,
					G4double z11, G4double z21, G4double z12, G4double z22)
{
  return (z11*(x2-xv)*(y2-yv) + z22*(xv-x1)*(yv-y1) + 
	  0.5*(z12*((x2-xv)*(yv-y1)+(xv-x1)*(y2-yv))+z21*((xv-x1)*(y2-yv)+(yv-y1)*(x2-xv))))
         / ((x2-x1)*(y2-y1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
