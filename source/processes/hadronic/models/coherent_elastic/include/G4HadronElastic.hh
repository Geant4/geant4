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
// $Id: G4HadronElastic.hh,v 1.9 2006-05-30 08:52:23 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Model: Low energy elastic scattering with 4-momentum balance 
// Derived fron G4LElastic of F.W. Jones, TRIUMF, 04-JUN-96
// Uses  G4ElasticHadrNucleusHE and G4VQCrossSection
//  
// Modified:
// 14-Dec-05 V.Ivanchenko rename the class
// 13-Apr-06 V.Ivanchenko move to coherent_elastic 
//
//
// Class Description
// Final state production model for hadron nuclear elastic scattering; 
// Class Description - End


#ifndef G4HadronElastic_h
#define G4HadronElastic_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

enum G4ElasticGenerator
{
  fLElastic = 0,
  fHElastic,
  fQElastic,
  fSWave
};

class G4ParticleDefinition;
class G4VQCrossSection;
class G4ElasticHadrNucleusHE;

class G4HadronElastic : public G4HadronicInteraction
{
public:

  G4HadronElastic(G4double elim = 100.*keV, 
		  G4double plow = 20.*MeV, 
		  G4double ehigh= GeV);

  virtual ~G4HadronElastic();
 
  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
				  G4Nucleus & targetNucleus);

  G4VQCrossSection* GetCS();

  G4ElasticHadrNucleusHE* GetHElastic();

  void SetMomentumLow(G4double value);

  void SetKinEnergyHigh(G4double value);

  G4double SampleT(G4double p, G4double m1, G4double m2, G4double A);

private:

  G4int Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
	     G4int iend,
	     G4double aa, G4double bb, G4double cc, G4double dd, 
	     G4double rr);

  G4double Fctcos(G4double t, 
		  G4double aa, G4double bb, G4double cc, G4double dd, 
		  G4double rr);

  void Defs1(G4double p, G4double px, G4double py, G4double pz, 
	     G4double pxinc, G4double pyinc, G4double pzinc, 
	     G4double* pxnew, G4double* pynew, G4double* pznew);

  G4VQCrossSection*           qCManager;
  G4ElasticHadrNucleusHE*     hElastic;

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theDeuteron;
  const G4ParticleDefinition* theAlpha;

  G4double ekinlim;  // in MeV
  G4double plablow;  // in MeV/c
  G4double ekinhigh;  // in MeV/c
};

inline void G4HadronElastic::SetMomentumLow(G4double value)
{
  plablow = value;
}

inline void G4HadronElastic::SetKinEnergyHigh(G4double value)
{
  ekinhigh = value;
}

#endif
