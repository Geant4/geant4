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
// GEANT4 physics class: G4QPhotoNuclearCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//

#ifndef G4QPhotoNuclearCrossSection_h
#define G4QPhotoNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4QVInelasticCrossSection.hh"
#include <vector>

class G4QPhotoNuclearCrossSection : public G4QVInelasticCrossSection
{
public:

  G4QPhotoNuclearCrossSection() : G4QVInelasticCrossSection() {}
  ~G4QPhotoNuclearCrossSection() {}

  G4double GetCrossSection(G4double Energy, G4int Z, G4int N);
  G4double ThresholdEnergy(G4int Z, G4int N);
  G4int    GetFunctions(G4double A, G4double* y, G4double* z);// y&z are pointers to arrays

// Body
private:
  static G4int     lastN;   // The last N of calculated nucleus
  static G4int     lastZ;   // The last Z of calculated nucleus
  static G4double  lastSig; // Last value of the Cross Section
  static G4double* lastGDR; // Pointer to the last array of GDR cross sections
  static G4double* lastHEN; // Pointer to the last array of HEn cross sections
  static G4double  lastE;   // Last used in the cross section Energy
  static G4double  lastTH;  // Last value of the Energy Threshold (A-dependent)
  static G4double  lastSP;  // Last value of the ShadowingPomeron (A-dependent)
};

// Gives the threshold energy for different nuclei (min of p- and n-threshold)
inline G4double G4QPhotoNuclearCrossSection::ThresholdEnergy(G4int Z, G4int N)
{
  static const G4double mNeut = G4NucleiProperties::GetNuclearMass(1,0);
  static const G4double mProt = G4NucleiProperties::GetNuclearMass(1,1);
  static const G4double mAlph = G4NucleiProperties::GetNuclearMass(4,2);
  // ---------
  static const G4double infEn = 9.e27;

  G4int A=Z+N;
  if(A<1) return infEn;
  else if(A==1) return 134.9766; // Pi0 threshold for the nucleon
  G4double mT = 0.;
  if(G4NucleiPropertiesTable::IsInTable(Z,A)) mT=G4NucleiProperties::GetNuclearMass(A,Z);
  else return 0.; // If it's not in the Table of Stable Nuclei, then unstable & Threshold=0
  // ---------
  G4double mP = infEn;
  if(Z && G4NucleiPropertiesTable::IsInTable(Z-1,A-1))
    mP=G4NucleiProperties::GetNuclearMass(A-1,Z-1); // Residual mass for a proton

  G4double mN = infEn;
  if(N && G4NucleiPropertiesTable::IsInTable(Z,A-1))
    mN=G4NucleiProperties::GetNuclearMass(A-1,Z);   // Residual mass for a neutron

  G4double mA= infEn;
  if(N>1&&Z>1&&G4NucleiPropertiesTable::IsInTable(Z-2,A-4))
	mA=G4NucleiProperties::GetNuclearMass(A-4,Z-2); // Residual mass for an alpha

  G4double dP = mP+mProt-mT, dN = mN+mNeut-mT, dA = mA + mAlph - mT;
  if(dP<dN) dN=dP; if(dA<dN) dN=dA; return dN;
}

#endif
