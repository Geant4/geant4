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
// $Id: G4PhotoNuclearCrossSection.hh,v 1.6 2002-05-22 11:48:45 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4PhotoNuclearCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//

#ifndef G4PhotoNuclearCrossSection_h
#define G4PhotoNuclearCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
/////////#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
//#include "G4QPDGCode.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "g4std/vector"

class G4PhotoNuclearCrossSection : public G4VCrossSectionDataSet
{
public:

  G4PhotoNuclearCrossSection()               // Constructor @@??
  {
	 //theHadronCrossSections = G4HadronCrossSections::Instance();
  }

  ~G4PhotoNuclearCrossSection() {}

  G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Element* anElement)
  {
	//return theHadronCrossSections->IsApplicable(aParticle, anElement);
	// Possible prototype
	G4bool result = false;
	if( aParticle->GetDefinition()->GetPDGEncoding()==22) result = true;
	return result;
  }

  G4double GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement,
                           G4double T=0.);
  //{
  //  return theHadronCrossSections->GetInelasticCrossSection(aParticle,
  //                                                          anElement);
  //}

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4int    GetFunctions(G4double a, G4double* y, G4double* z);
  //G4double LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN);
  G4double EquLinearFit(G4double X, G4int N,const G4double X0,const G4double XD, const G4double* Y);
  G4double ThresholdEnergy(G4int Z, G4int N);

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

  //G4HadronCrossSections* theHadronCrossSections;
};

// Gives the threshold energy for different nuclei (min of p- and n-threshold)
inline G4double G4PhotoNuclearCrossSection::ThresholdEnergy(G4int Z, G4int N)
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
  else
  {
    G4cerr<<"G4PhotoNucCrossSect.hh::ThreshEn:Z="<<Z<<",A="<<A<<" element isn't in G4NucPr"<<G4endl;
    return 0.;                // If it is not in the Table of Stable Nuclei, then the Threshold=inf
  }
  // ---------
  G4double mP= infEn;
  //if(Z) mP= G4QPDGCode(111).GetNuclMass(Z-1,N,0);
  if(Z&&G4NucleiPropertiesTable::IsInTable(Z-1,A-1)) mP=G4NucleiProperties::GetNuclearMass(A-1,Z-1);
  else
  {
    G4cerr<<"G4PhotoNucCrossSect.hh::ThrEn:Z="<<Z-1<<",A="<<A-1<<" element isn't in G4NucP"<<G4endl;
  }
  G4double mN= infEn;
  //if(N) mN= G4QPDGCode(111).GetNuclMass(Z,N-1,0);
  if(N&&G4NucleiPropertiesTable::IsInTable(Z,A-1)) mN=G4NucleiProperties::GetNuclearMass(A-1,Z);
  else
  {
    G4cerr<<"G4PhotoNucCrossSect.hh::ThreshEn:Z="<<Z<<",A="<<A-1<<" element isn't in G4NuP"<<G4endl;
  }
  G4double dP= mP+mProt-mT;
  G4double dN= mN+mNeut-mT;
  if(dP<dN)dN=dP;
  return dN;
}

#endif
