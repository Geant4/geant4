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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ElectroNuclearCrossSection.hh,v 1.6 2001-11-30 14:58:04 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectroNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
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
						   G4double temperature=0.);

  G4double GetEffectivePhotonEnergy();

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4int    GetFunctions(G4double a, G4double* y, G4double* z);
  G4double LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN);
  G4double ThresholdEnergy(G4int Z, G4int N);
  G4double HighEnergyPhi(G4double lE);
  G4double HighEnergyFun(G4double lE);
  G4double SolveTheEquation(G4double f);
  G4double Fun(G4double x);
  G4double DFun(G4double x);

// Body
private:
  static G4int     lastF;    // Last used in the cross section TheFirstBin
  static G4double* lastPhi;  // Pointer to the last array of the Phi function
  static G4double* lastFun;  // Pointer to the last array of the Fun function
  static G4int     lastL;    // Last used in the cross section TheLastBin
  static G4double  lastLE;   // Last used in the cross section TheLogE
  static G4double  lastSig;  // Last value of the Cross Section
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

inline G4double G4ElectroNuclearCrossSection::DFun(G4double x)
{
  static const G4double f21=.75;
  static const G4double f22=.78;
  static const G4double f23=.9;
  G4double f2=1.;
  if     (lastH<.005) f2=f21;
  else if(lastH<0.01) f2=f22;
  else if(lastH<0.07) f2=f23;
  return (lastLE-x)*(0.0116*exp(0.16*x)+f2*exp(-0.26*x));
}

inline G4double G4ElectroNuclearCrossSection::Fun(G4double x)
{return (lastLE*HighEnergyPhi(x)-HighEnergyFun(x));}

inline G4double G4ElectroNuclearCrossSection::HighEnergyPhi(G4double lEn)
{
  static const G4double le=log(2000.);
  static const G4double c1=0.16;
  static const G4double f1=.0116/c1;
  static const G4double e1=exp(c1*le);
  static const G4double c2=-0.26;
  static const G4double f21=.75/c2;
  static const G4double f22=.78/c2;
  static const G4double f23=.9/c2;
  static const G4double f24=1./c2;
  static const G4double e2=exp(c2*le);
  G4double f2=f24;
  if     (lastH<.005) f2=f21;
  else if(lastH<0.01) f2=f22;
  else if(lastH<0.07) f2=f23;
  return f1*(exp(c1*lEn)-e1)+f2*(exp(c2*lEn)-e2);
}

inline G4double G4ElectroNuclearCrossSection::HighEnergyFun(G4double lEn)
{
  static const G4double le=log(2000.);
  static const G4double c1=0.16;
  static const G4double f1=.0116/c1/c1;
  static const G4double e1=(c1*le-1.)*exp(c1*le);
  static const G4double c2=-0.26;
  static const G4double f21=.75/c2/c2;
  static const G4double f22=.78/c2/c2;
  static const G4double f23=.9/c2/c2;
  static const G4double f24=1./c2/c2;
  static const G4double e2=(c2*le-1.)*exp(c2*le);
  G4double f2=f24;
  if     (lastH<.005) f2=f21;
  else if(lastH<0.01) f2=f22;
  else if(lastH<0.07) f2=f23;
  G4double h1=c1*lEn;
  G4double h2=c2*lEn;
  return f1*((h1-1.)*exp(h1)-e1)+f2*((h2-1.)*exp(h2)-e2);
}

#endif
