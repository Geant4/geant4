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
// $Id: G4PhotoNuclearCrossSection.hh,v 1.3 2001-10-26 13:12:25 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4PhotoNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
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
						   G4double temperature=0.);
  //{
  //  return theHadronCrossSections->GetInelasticCrossSection(aParticle,
  //                                                          anElement);
  //}

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4double GetGDRc1(G4int Z, G4int N);
  G4double GetGDRp1(G4int Z, G4int N);
  G4double GetGDRt1(G4int Z, G4int N);
  G4double GetGDRs1(G4int Z, G4int N);
  G4double GetGDRc2(G4int Z, G4int N);
  G4double GetGDRp2(G4int Z, G4int N);
  G4double GetGDRt2(G4int Z, G4int N);
  G4double GetGDRs2(G4int Z, G4int N);
  G4double GetQDAmp(G4int Z, G4int N);
  G4double GetDelAm(G4int Z, G4int N);
  G4double GetDelWd(G4int Z, G4int N);
  G4double GetDelPs(G4int Z, G4int N);
  G4double GetDelTh(G4int Z, G4int N);
  G4double GetDelSl(G4int Z, G4int N);
  G4double GetRopAm(G4int Z, G4int N);
  G4double GetRopWd(G4int Z, G4int N);
  G4double GetRopPs(G4int Z, G4int N);
  G4double LinearFit(G4double X, G4int N, const G4double* XN, const G4double* YN);
  G4double ThresholdEnergy(G4int Z, G4int N);

// Body
//private:

  //G4HadronCrossSections* theHadronCrossSections;
};

// Calculate the logAmplitude of the 1-st GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRc1(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={4.2,13.9,13.9,13.6,20.5,28.2,28.7,28.5,29.,28.4,28.15,27.8,25.9};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the A-power of the 1-st GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRp1(G4int Z, G4int N)
{
  G4double p=8.;
  G4int A=Z+N;
  if(A<12) p=6.;
  if(A< 8) p=4.;
  if(A< 4) p=2.;
  return p;
}

// Calculate the Threshold of the 1-st GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRt1(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={1.4,3.13,3.08,2.9,3.09,3.09,3.09,3.02,2.98,2.9,2.745,2.585,2.42};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the Slope of the 1-st GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRs1(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={.12,.12,.12,.12,.06,.03,.03,.06,.05,.065,.06,.059,.061};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the logAmplitude of the 2-nd GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRc2(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={1.85,7.5,6.3,8.2,12.35,15.8,16.1,16.2,16.8,17.1,16.1,15.5,16.6};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the A-power of the 2-nd GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRp2(G4int Z, G4int N)
{
  G4double p=4.;
  G4int A=Z+N;
  if(A<12) p=3.;
  if(A< 8) p=2.;
  if(A< 4) p=1.;
  return p;
}

// Calculate the Threshold of the 2-nd GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRt2(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={1.4,3.22,3.11,3.39,3.48,3.34,3.46,3.35,3.4,3.22,3.09,3.05,2.6};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the Slope of the 2-nd GDR maximum
inline G4double G4PhotoNuclearCrossSection::GetGDRs2(G4int Z, G4int N)
{
  static const G4int nN=13;
  static G4double X[nN]={0.693,1.386,1.792,1.946,2.197,2.485,2.773,3.296,3.689,4.152,4.777,5.334,
						 5.472};
  static G4double Y[nN]={.12,.094,.09,.088,.14,.082,.079,.074,.071,.065,.061,.058,.05};

  return LinearFit(log(Z+N), nN, X, Y);
}

// Calculate the Amplitude of the QuasiDeuteron region [exp/(1+exp)]
inline G4double G4PhotoNuclearCrossSection::GetQDAmp(G4int Z, G4int N)
{
  G4double A=Z+N;
  G4double lnA=log(A);
  return exp(-1.7+lnA*0.84)/(1.+exp(7*(2.38-lnA)));
}

// Calculate the Amplitude of the Delta Resonance [.41*(Z+N)]
inline G4double G4PhotoNuclearCrossSection::GetDelAm(G4int Z, G4int N)
{
  G4double A=Z+N;
  return .41*A;
}

// Calculate the Width of the Delta Resonance [11.9-ln(A)*1.24]
inline G4double G4PhotoNuclearCrossSection::GetDelWd(G4int Z, G4int N)
{
  G4double A=Z+N;
  G4double lnA=log(A);
  return 11.9-lnA*1.24;
}

// Calculate the Position of the Delta Resonance [5.84-.09/(1+.003*A*A)]
inline G4double G4PhotoNuclearCrossSection::GetDelPs(G4int Z, G4int N)
{
  G4double A=Z+N;
  return 5.84-.09/(1+.003*A*A);
}

// Calculate the Threshold of the Delta Resonance [5.13-.00075*A]
inline G4double G4PhotoNuclearCrossSection::GetDelTh(G4int Z, G4int N)
{
  G4double A=Z+N;
  return 5.13-0.00075*A;
}

// Calculate the Threshold of the Delta Resonance [.04->.09]
inline G4double G4PhotoNuclearCrossSection::GetDelSl(G4int Z, G4int N)
{
  G4double A=Z+N;
  if(A<7) return .04;
  return .09;
}

// Calculate the Amplitude of the Roper Resonance [-2.+ln(A)*0.84]
inline G4double G4PhotoNuclearCrossSection::GetRopAm(G4int Z, G4int N)
{
  G4double A=Z+N;
  G4double lnA=log(A);
  return exp(-2.+lnA*0.84);
}

// Calculate the Width of the Roper Resonance [.1+1.65*ln(A)]
inline G4double G4PhotoNuclearCrossSection::GetRopWd(G4int Z, G4int N)
{
  G4double A=Z+N;
  G4double lnA=log(A);
  return .1+1.65*lnA;
}

// Calculate the Position of the Roper Resonance [6.46+.061*ln(A)]
inline G4double G4PhotoNuclearCrossSection::GetRopPs(G4int Z, G4int N)
{
  G4double A=Z+N;
  G4double lnA=log(A);
  return 6.46+.061*lnA;
}


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

#endif
