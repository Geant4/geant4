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
//  Calculation of the total, elastic and inelastic cross-sections
//  of anti-nucleon and anti-nucleus interactions with nuclei 
//  based on Glauber approach,  V. Grishine formulaes for 
//  interpolations (ref. V.M.Grichine, Eur.Phys.J., C62(2009) 399;   
//  NIM, B267 (2009) 2460) and our parametrization of hadron-nucleon 
//  cross-sections 
// 
//
//   Created by A.Galoyan and V. Uzhinsky, 18.11.2010  


#ifndef G4ComponentAntiNuclNuclearXS_h
#define G4ComponentAntiNuclNuclearXS_h

#include <CLHEP/Units/PhysicalConstants.h>  // pi, fermi,..

#include "globals.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiHe3.hh"
#include "G4AntiTriton.hh"
#include "G4AntiAlpha.hh"
#include "G4Nucleus.hh"
#include "G4Pow.hh"
#include "G4VComponentCrossSection.hh"


class G4ParticleDefinition;


class G4ComponentAntiNuclNuclearXS : public G4VComponentCrossSection {

  public:
    G4ComponentAntiNuclNuclearXS ();
    virtual ~G4ComponentAntiNuclNuclearXS ();
    virtual G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				                 G4double kinEnergy, G4int Z, G4int A);
    virtual G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				                 G4double kinEnergy, G4int Z, G4double A);
    virtual G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					             G4double kinEnergy, G4int Z, G4int A);
    virtual G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					             G4double kinEnergy, G4int Z, G4double A);
    virtual G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					           G4double kinEnergy, G4int Z, G4double A);
    virtual G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					           G4double kinEnergy, G4int Z, G4int A);
    virtual void BuildPhysicsTable(const G4ParticleDefinition&) {}
    virtual void DumpPhysicsTable(const G4ParticleDefinition&)  {}
    virtual void CrossSectionDescription(std::ostream&) const;
    // Method for calculation of Anti-Hadron Nucleon Total Cross-section
    G4double GetAntiHadronNucleonTotCrSc(const G4ParticleDefinition* aParticle, G4double kinEnergy);
    // Method for calculation of Anti-Hadron Nucleon Elastic Cross-section
    G4double GetAntiHadronNucleonElCrSc(const G4ParticleDefinition* aParticle, G4double kinEnergy);

  private:
    G4double fRadiusEff;  // Effective Radius for AntiNucleus 
    G4double fTotalXsc, fElasticXsc, fInelasticXsc;
    G4double fAntiHadronNucleonTotXsc, fAntiHadronNucleonElXsc; 
    G4double Elab, S, SqrtS ;
    G4double Mn, b0, b2,  SqrtS0, S0, R0;  // Parameters for AntiHadron-Nucleon Xsc  
    G4ParticleDefinition* theAProton;
    G4ParticleDefinition* theANeutron;
    G4ParticleDefinition* theADeuteron;
    G4ParticleDefinition* theATriton;
    G4ParticleDefinition* theAAlpha;
    G4ParticleDefinition* theAHe3;

    const G4double ReffTot[5][5] =  { {0.000, 3.800, 3.300, 3.300, 2.376},    // Pbar   + p, d, t, He3, He4
                                      {3.800, 3.238, 3.144, 3.144, 2.544},    // Dbar   + p, d, t, He3, He4
                                      {3.300, 3.144, 3.075, 3.075, 3.589},    // Tbar   + p, d, t, He3, He4
                                      {3.300, 3.144, 3.075, 3.075, 2.589},    // He3bar + p, d, t, He3, He4
                                      {2.376, 2.544, 3.589, 3.598, 2.241} };  // He4bar + p, d, t, He3, He4
    const G4double ReffInel[5][5] = { {0.000, 3.582, 3.105, 3.105, 2.209},    // Pbar   + p, d, t, He3, He4
                                      {3.582, 3.169, 3.066, 3.066, 2.498},    // Dbar   + p, d, t, He3, He4
                                      {3.105, 3.066, 2.973, 2.973, 2.508},    // Tbar   + p, d, t, He3, He4
                                      {3.105, 3.066, 2.973, 2.973, 2.508},    // He3bar + p, d, t, He3, He4
                                      {2.209, 2.498, 2.508, 2.508, 2.158} };  // He4bar + p, d, t, He3, He4

    const G4Pow* theG4Pow = G4Pow::GetInstance();
};

#endif
