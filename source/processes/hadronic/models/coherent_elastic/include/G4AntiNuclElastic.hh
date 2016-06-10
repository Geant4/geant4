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
// $Id: G4AntiNuclElastic.hh, 2011-04-12 agaloyan
//
// Geant4 Header : G4AntiNuclElastic
//
// Author : A. Galoyan  
//  
// Class Description
// Model for AntiNuclear Nuclear Elastic Scattering;
// Class Description - End

#ifndef G4AntiNuclElastic_h
#define G4AntiNuclElastic_h 1
 
#include "G4HadronElastic.hh"
#include "globals.hh"
#include "G4Nucleus.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
	
class G4ParticleDefinition;

class G4AntiNuclElastic : public G4HadronElastic
{
public:

  G4AntiNuclElastic();

  virtual ~G4AntiNuclElastic();
 
  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);

 G4double SampleThetaCMS(const G4ParticleDefinition* p, G4double plab, 
 G4int Z, G4int A);

 G4double SampleThetaLab(const G4ParticleDefinition* p, 
                                G4double plab, G4int Z, G4int A);

 
 G4double CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    );
 
 G4double CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 );

 G4double CalculateAm( G4double momentum, G4double n, G4double Z);
 
 G4double DampFactor(G4double z);

 G4double BesselJzero(G4double z);

 G4double BesselJone(G4double z);

 G4double BesselOneByArg(G4double z); 

 G4double GetcosTeta1( G4double plab, G4int A);

 inline G4ComponentAntiNuclNuclearXS* GetComponentCrossSection();
 
private:

  // Assignment operator and copy constructor
  G4AntiNuclElastic & operator=(const G4AntiNuclElastic &right);
  G4AntiNuclElastic(const G4AntiNuclElastic&);

  G4ComponentAntiNuclNuclearXS* cs;    //cross section of antiA-A interaction

  G4double fTetaCMS;        //  sampled Theta in CMS 
  G4double fThetaLab;        //sampled Theta in Lab system
  const G4ParticleDefinition*  fParticle;
  G4double fWaveVector;
  G4double fBeta;         // velosity of projectile 
  G4double fZommerfeld;   // parameter of Zommerfeld for calculation of Coulomb cross-section
  G4double fAm;           // parameter for calculation of Coulomb cross-section
  G4double fRa;           // Radius of target
  G4double fRef;          // Effective radiuse for Calculation of hadron cross-section
  G4double fceff;         //  Effective diffuse parameter

  G4ThreeVector fbst;          // boost vector
  G4double fptot;         // momentum of projectile in CMS system
  G4double fTmax;          

  G4ParticleDefinition* theAProton;
  G4ParticleDefinition* theANeutron;
  G4ParticleDefinition* theADeuteron;
  G4ParticleDefinition* theATriton;
  G4ParticleDefinition* theAAlpha;
  G4ParticleDefinition* theAHe3;
 
  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theDeuteron;
  G4ParticleDefinition* theAlpha;
    
};

inline G4ComponentAntiNuclNuclearXS* 
G4AntiNuclElastic::GetComponentCrossSection()
{
  return cs;
}

#endif


