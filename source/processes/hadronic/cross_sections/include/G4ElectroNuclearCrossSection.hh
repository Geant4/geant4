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
// $Id: G4ElectroNuclearCrossSection.hh,v 1.8 2002-12-12 19:16:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectroNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//

#ifndef G4ElectroNuclearCrossSection_h
#define G4ElectroNuclearCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "g4std/vector"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4ElectroNuclearCrossSection : public G4VCrossSectionDataSet
{
public:

  G4ElectroNuclearCrossSection()               // Constructor @@??
  {
	 //theHadronCrossSections = G4HadronCrossSections::Instance();
  }

  ~G4ElectroNuclearCrossSection() {}

  G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Element* anElement)
  {
	//return theHadronCrossSections->IsApplicable(aParticle, anElement);
	// Possible prototype
	G4bool result = false;
	if( aParticle->GetDefinition()==G4Electron::ElectronDefinition()) result = true;
	if( aParticle->GetDefinition()==G4Positron::PositronDefinition()) result = true;
	return result;
  }

  G4double GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement,
                           G4double T=0.);

  G4double GetEquivalentPhotonEnergy();

  G4double GetVirtualFactor(G4double nu, G4double Q2);

  G4double GetEquivalentPhotonQ2(G4double nu);

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4int    GetFunctions(G4double a, G4double* x, G4double* y, G4double* z);
  //G4double LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN);
  G4double ThresholdEnergy(G4int Z, G4int N);
  G4double HighEnergyJ1(G4double lE);
  G4double HighEnergyJ2(G4double lE);
  G4double HighEnergyJ3(G4double lE);
  G4double SolveTheEquation(G4double f);
  G4double Fun(G4double x);
  G4double DFun(G4double x);

// Body
private:
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4int     lastF;    // Last used in the cross section TheFirstBin
  static G4double* lastJ1;   // Pointer to the last array of the J1 function
  static G4double* lastJ2;   // Pointer to the last array of the J2 function
  static G4double* lastJ3;   // Pointer to the last array of the J3 function
  static G4int     lastL;    // Last used in the cross section TheLastBin
  static G4double  lastE;    // Last used in the cross section Energy
  static G4double  lastTH;   // Last value of the Energy Threshold
  static G4double  lastSig;  // Last value of the Cross Section
  static G4double  lastG;    // Last value of gamma=lnE-ln(me)
  static G4double  lastH;    // Last value of the High energy A-dependence
};

// Gives the threshold energy for different nuclei (min of p- and n-threshold)
inline G4double G4ElectroNuclearCrossSection::ThresholdEnergy(G4int Z, G4int N)
{
  // CHIPS - Direct GEANT
  //static const G4double mNeut = G4QPDGCode(2112).GetMass();
  //static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mNeut = G4NucleiProperties::GetNuclearMass(1,0);
  static const G4double mProt = G4NucleiProperties::GetNuclearMass(1,1);
  // ---------
  static const G4double infEn = 9.e27;

  G4int A=Z+N;
  if(A<1) return infEn;
  else if(A==1) return 134.9766; // Pi0 threshold for the nucleon
  // CHIPS - Direct GEANT
  //G4double mT= G4QPDGCode(111).GetNuclMass(Z,N,0);
  G4double mT= 0.;
  if(G4NucleiPropertiesTable::IsInTable(Z,A)) mT=G4NucleiProperties::GetNuclearMass(A,Z);
  else return 0.;                // If it is not in the Table of Stable Nuclei, then the Threshold=0
  // ---------
  G4double mP= infEn;
  //if(Z) mP= G4QPDGCode(111).GetNuclMass(Z-1,N,0);
  if(Z&&G4NucleiPropertiesTable::IsInTable(Z-1,A-1)) mP=G4NucleiProperties::GetNuclearMass(A-1,Z-1);
  else return infEn;
  G4double mN= infEn;
  //if(N) mN= G4QPDGCode(111).GetNuclMass(Z,N-1,0);
  if(N&&G4NucleiPropertiesTable::IsInTable(Z,A-1)) mN=G4NucleiProperties::GetNuclearMass(A-1,Z);
  else return infEn;
  G4double dP= mP+mProt-mT;
  G4double dN= mN+mNeut-mT;
  if(dP<dN)dN=dP;
  return dN;
}

inline G4double G4ElectroNuclearCrossSection::DFun(G4double x) // Original PhoNuc cross section
{
  static const G4double shd=1.0734;                    // HE PomShadowing(D)
  static const G4double poc=0.0375;                    // HE Pomeron coefficient
  static const G4double pos=16.5;                      // HE Pomeron shift
  static const G4double reg=.11;                       // HE Reggeon slope
  static const G4double mel=0.5109989;                 // Mass of electron in MeV
  static const G4double lmel=log(mel);                 // Log of electron mass
  G4double lE=lastG+lmel;
  return poc*(lE-pos)+shd*exp(-reg*lE);
}

inline G4double G4ElectroNuclearCrossSection::Fun(G4double x) // Integrated PhoNuc cross section
{
  static const G4double mel=0.5109989;                 // Mass of electron in MeV
  static const G4double lmel=log(mel);                 // Log of electron mass
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  G4double lE=lastG+lmel;
  G4double HE2=HighEnergyJ2(lE);
  return dlg1*HighEnergyJ1(lE)-lgoe*(HE2+HE2-HighEnergyJ3(lE)/lastE);
}

inline G4double G4ElectroNuclearCrossSection::HighEnergyJ1(G4double lEn)
{
  static const G4double le=log(50000.); // log(E0)
  static const G4double le2=le*le;      // log(E0)^2
  static const G4double a=.0375;        // a
  static const G4double ha=a*.5;        // a/2
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=0.11;         // d
  static const G4double cd=1.0734/d;    // c/d
  static const G4double ele=exp(-d*le); // E0^(-d)
  return ha*(lEn*lEn-le2)-ab*(lEn-le)-cd*(exp(-d*lEn)-ele);
}

inline G4double G4ElectroNuclearCrossSection::HighEnergyJ2(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=log(e);      // log(E0)
  static const G4double le1=(le-1.)*e;  // (log(E0)-1)*E0
  static const G4double a=.0375;        // a
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=1.-0.11;      // 1-d
  static const G4double cd=1.0734/d;    // c/(1-d)
  static const G4double ele=exp(d*le);  // E0^(1-d)
  return a*((lEn-1.)*lastE-le1)-ab*(lastE-e)+cd*(exp(d*lEn)-ele);
}

inline G4double G4ElectroNuclearCrossSection::HighEnergyJ3(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=log(e);      // log(E0)
  static const G4double e2=e*e;         // E0^2
  static const G4double leh=(le-.5)*e2; // (log(E0)-.5)*E0^2
  static const G4double ha=.0375*.5;    // a/2
  static const G4double hab=ha*16.5;    // a*b/2
  static const G4double d=2.-.11;       // 2-d
  static const G4double cd=1.0734/d;    // c/(2-d)
  static const G4double ele=exp(d*le);  // E0^(2-d)
  G4double lastE2=lastE*lastE;
  return ha*((lEn-.5)*lastE2-leh)-hab*(lastE2-e2)+cd*(exp(d*lEn)-ele);
}

#endif
