//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4HadronElastic.hh,v 1.32 2010-01-13 15:42:06 mkossov Exp $
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
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy, below which S-wave is sampled
// 20-Oct-06 V.Ivanchenko default ekinhigh = GeV (use HE model)
// 16-Nov-06 V.Ivanchenko remove definition 
// 16-Nov-06 V.Ivanchenko default ekinhigh = 0.4*GeV 
// 16-Nov-06 V.Ivanchenko cleanup and rename Set methods and variables 
// 28-Mar-07 V.Ivanchenko add NIST manager
// 11-May-07 V.Ivanchenko remove unused method Defs1
// 13.01.10: M.Kosov: Use G4Q(Pr/Neut)ElasticCS instead of G4QElasticCS
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

  G4HadronElastic(G4ElasticHadrNucleusHE* HModel = 0);

  virtual ~G4HadronElastic();
 
  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
				  G4Nucleus & targetNucleus);

  G4VQCrossSection* GetCS();

  G4ElasticHadrNucleusHE* GetHElastic();

  void SetPlabLowLimit(G4double value);

  void SetHEModelLowLimit(G4double value);

  void SetQModelLowLimit(G4double value);

  void SetLowestEnergyLimit(G4double value);

  void SetRecoilKinEnergyLimit(G4double value);

  G4double SampleT(G4double p, G4double m1, G4double m2, G4double A);

private:

  G4int Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
	     G4int iend,
	     G4double aa, G4double bb, G4double cc, G4double dd, 
	     G4double rr);

  G4double Fctcos(G4double t, 
		  G4double aa, G4double bb, G4double cc, G4double dd, 
		  G4double rr);

  static G4VQCrossSection* pCManager;
  static G4VQCrossSection* nCManager;

  G4ElasticHadrNucleusHE*     hElastic;

  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theDeuteron;
  G4ParticleDefinition* theAlpha;
  const G4ParticleDefinition* thePionPlus;
  const G4ParticleDefinition* thePionMinus;

  G4double lowEnergyRecoilLimit;  
  G4double lowEnergyLimitHE;  
  G4double lowEnergyLimitQ;  
  G4double lowestEnergyLimit;  
  G4double plabLowLimit;

};

inline void G4HadronElastic::SetRecoilKinEnergyLimit(G4double value)
{
  lowEnergyRecoilLimit = value;
}

inline void G4HadronElastic::SetPlabLowLimit(G4double value)
{
  plabLowLimit = value;
}

inline void G4HadronElastic::SetHEModelLowLimit(G4double value)
{
  lowEnergyLimitHE = value;
}

inline void G4HadronElastic::SetQModelLowLimit(G4double value)
{
  lowEnergyLimitQ = value;
}

inline void G4HadronElastic::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}

#endif
