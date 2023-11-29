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
//
#ifndef G4FTFParameters_h
#define G4FTFParameters_h 1

#include <CLHEP/Units/SystemOfUnits.h>
#include <vector>
#include "G4Types.hh"
#include "G4Exp.hh"
#include "G4FTFTunings.hh"

class G4ParticleDefinition;
class G4VComponentCrossSection;
class G4LundStringFragmentation;


class G4FTFParameters {
  public:
    G4FTFParameters();
    ~G4FTFParameters();

    void InitForInteraction( const G4ParticleDefinition* , G4int theA, G4int theZ, G4double s );

    // Set geometrical parameteres
    void SethNcmsEnergy( const G4double s );
    void SetTotalCrossSection( const G4double Xtotal );
    void SetElastisCrossSection( const G4double Xelastic );
    void SetInelasticCrossSection( const G4double Xinelastic );
    void SetProbabilityOfElasticScatt( const G4double Xtotal, const G4double Xelastic );
    void SetProbabilityOfElasticScatt( const G4double aValue );
    void SetProbabilityOfAnnihilation( const G4double aValue ); 
    void SetRadiusOfHNinteractions2( const G4double Radius2 );

    void SetSlope( const G4double Slope );
    void SetGamma0( const G4double Gamma0 );
    G4double GammaElastic( const G4double impactsquare );

    // Set parameters of elastic scattering
    void SetAvaragePt2ofElasticScattering( const G4double aPt2 );

    // Set parameters of excitations
    void SetParams( const G4int ProcN, 
                    const G4double A1, const G4double B1, const G4double A2, const G4double B2,
                    const G4double A3, const G4double Atop, const G4double Ymin );

    void SetDeltaProbAtQuarkExchange( const G4double aValue );
    void SetProbOfSameQuarkExchange( const G4double aValue );

    void SetProjMinDiffMass( const G4double aValue );
    void SetProjMinNonDiffMass( const G4double aValue );
    //void SetProbabilityOfProjDiff( const G4double aValue );
    void SetProbLogDistrPrD( const G4double aValue );

    void SetTarMinDiffMass( const G4double aValue ); 
    void SetTarMinNonDiffMass( const G4double aValue );
    //void SetProbabilityOfTarDiff( const G4double aValue );

    void SetAveragePt2( const G4double aValue );
    void SetProbLogDistr( const G4double aValue );

    // Set parameters of a string kink
    void SetPt2Kink( const G4double aValue );
    void SetQuarkProbabilitiesAtGluonSplitUp( const G4double Puubar, const G4double Pddbar, 
                                              const G4double Pssbar );

    // Set parameters of nuclear destruction
    void SetMaxNumberOfCollisions( const G4double aValue, const G4double bValue );
    void SetProbOfInteraction( const G4double aValue );

    void SetCofNuclearDestructionPr( const G4double aValue );
    void SetCofNuclearDestruction( const G4double aValue );
    void SetR2ofNuclearDestruction( const G4double aValue );

    void SetExcitationEnergyPerWoundedNucleon( const G4double aValue );

    void SetDofNuclearDestruction( const G4double aValue );
    void SetPt2ofNuclearDestruction( const G4double aValue );
    void SetMaxPt2ofNuclearDestruction( const G4double aValue );

    // Get geometrical parameteres
    G4double GetTotalCrossSection();
    G4double GetElasticCrossSection();
    G4double GetInelasticCrossSection();

    G4double GetProbabilityOfInteraction( const G4double impactsquare );
    G4double GetInelasticProbability( const G4double impactsquare );
    G4double GetProbabilityOfElasticScatt();
    G4double GetSlope();
    G4double GetProbabilityOfAnnihilation(); 

    // Get parameters of elastic scattering
    G4double GetAvaragePt2ofElasticScattering();

    // Get parameters of excitations
    G4double GetProcProb( const G4int ProcN, const G4double y );

    G4double GetDeltaProbAtQuarkExchange();
    G4double GetProbOfSameQuarkExchange();

    G4double GetProjMinDiffMass();
    G4double GetProjMinNonDiffMass();
    G4double GetProbLogDistrPrD();

    G4double GetTarMinDiffMass();
    G4double GetTarMinNonDiffMass();

    G4double GetAveragePt2();
    G4double GetProbLogDistr();

    // Get parameters of a string kink
    G4double GetPt2Kink();
    std::vector< G4double > GetQuarkProbabilitiesAtGluonSplitUp();

    // Get parameters of nuclear destruction
    G4double GetMaxNumberOfCollisions();
    G4double GetProbOfInteraction();

    G4double GetCofNuclearDestructionPr();
    G4double GetCofNuclearDestruction();
    G4double GetR2ofNuclearDestruction();

    G4double GetExcitationEnergyPerWoundedNucleon();

    G4double GetDofNuclearDestruction();
    G4double GetPt2ofNuclearDestruction();
    G4double GetMaxPt2ofNuclearDestruction();

    // JVY, July 31, 2017: Is there any reason for NOT making 
    //                     all the members data private ???
    //
    // private: 

    // Initial energy of hN interactions
    G4double FTFhNcmsEnergy;  // Initial hN CMS energy

    // Geometrical parameteres
    G4double FTFXtotal;                      // Total X in mb
    G4double FTFXelastic;                    // Elastic X in mb
    G4double FTFXinelastic;                  // Inelastic X in mb
    G4double FTFXannihilation;               // Annihilation X in mb 
    G4double ProbabilityOfAnnihilation;      // Xannih/Xinelast     
    G4double ProbabilityOfElasticScatt;      // Xel/Xtot
    G4double RadiusOfHNinteractions2;        // Xtot/pi, in fm^2
    G4double FTFSlope;                       // in fm^-1
    G4double AvaragePt2ofElasticScattering;  // in MeV^2
    G4double FTFGamma0;

    // Parameters of excitations
    G4double ProcParams[5][7];

    G4double DeltaProbAtQuarkExchange;
    G4double ProbOfSameQuarkExchange;

    G4double ProjMinDiffMass;
    G4double ProjMinNonDiffMass;
    G4double ProbLogDistrPrD;
    G4double TarMinDiffMass;  
    G4double TarMinNonDiffMass;

    G4double AveragePt2;
    G4double ProbLogDistr;

    // Parameters of kink
    G4double Pt2kink;
    std::vector< G4double > QuarkProbabilitiesAtGluonSplitUp;

    // Parameters of nuclear destruction
    G4double MaxNumberOfCollisions;
    G4double ProbOfInelInteraction;

    G4double CofNuclearDestructionPr; // Cnd of nuclear destruction of projectile nucleus
    G4double CofNuclearDestruction;   // Cnd of nuclear destruction
    G4double R2ofNuclearDestruction;  // R2nd

    G4double ExcitationEnergyPerWoundedNucleon;

    G4double DofNuclearDestruction;       // Dispersion for momentum sampling
    G4double Pt2ofNuclearDestruction;     // Pt2
    G4double MaxPt2ofNuclearDestruction;  // Max Pt2

    G4bool EnableDiffDissociationForBGreater10; ///< Control over whether to do nucleon-hadron diffractive dissociation or not.

  private:
    G4LundStringFragmentation* StringMass;
    G4double GetMinMass( const G4ParticleDefinition* aParticle );

    void Reset(); 

    // Different sets of parameters (called "tunes") of the FTF model are possible.
    // These tunes are kept as std::array - instead of std::vector - members of this class,
    // because their size is fixed during a run, and expected to be small.
    // For the time being, separate parameters are kept for "baryons", "pions", and
    // the rest of "mesons"; if in the future we make more distinctions between
    // projectile types (e.g. kaons, anti-baryon, hyperons, etc.), then corresponding
    // new arrays will be introduced. In all cases, the size of these arrays is the
    // same (and kept as a static constant in the singleton G4FTFTunings).
    std::array< G4FTFParamCollBaryonProj, G4FTFTunings::sNumberOfTunes > fArrayParCollBaryonProj;
    std::array< G4FTFParamCollMesonProj,  G4FTFTunings::sNumberOfTunes > fArrayParCollMesonProj;
    std::array< G4FTFParamCollPionProj,   G4FTFTunings::sNumberOfTunes > fArrayParCollPionProj;

    // Glauber-Gribov hN x-section
    G4VComponentCrossSection* csGGinstance;
};


inline G4double G4FTFParameters::GammaElastic( const G4double impactsquare ) {
  return ( FTFGamma0 * G4Exp( -FTFSlope * impactsquare ) );
}

inline void G4FTFParameters::SethNcmsEnergy( const G4double S ) { 
  FTFhNcmsEnergy = S;
}

// Set geometrical parameteres

inline void G4FTFParameters::SetTotalCrossSection( const G4double Xtotal ) {
  FTFXtotal = Xtotal;
}

inline void G4FTFParameters::SetElastisCrossSection( const G4double Xelastic ) {
  FTFXelastic = Xelastic;
}

inline void G4FTFParameters::SetInelasticCrossSection( const G4double Xinelastic ) {
  FTFXinelastic = Xinelastic;
}

inline void G4FTFParameters::SetProbabilityOfElasticScatt( const G4double Xtotal, 
                                                           const G4double Xelastic ) { 
  if ( Xtotal == 0.0 ) {
    ProbabilityOfElasticScatt = 0.0;
  } else {
    ProbabilityOfElasticScatt = Xelastic / Xtotal;
  }
} 

inline void G4FTFParameters::SetProbabilityOfElasticScatt( const G4double aValue ) {
  ProbabilityOfElasticScatt = aValue;
}

inline void G4FTFParameters::SetProbabilityOfAnnihilation( const G4double aValue ) {
  ProbabilityOfAnnihilation = aValue;
}

inline void G4FTFParameters::SetRadiusOfHNinteractions2( const G4double Radius2 ) {
  RadiusOfHNinteractions2 = Radius2;
}

inline void G4FTFParameters::SetSlope( const G4double Slope ) {
  FTFSlope = 12.84 / Slope; // Slope is in GeV^-2, FTFSlope in fm^-2
} 

inline void G4FTFParameters::SetGamma0( const G4double Gamma0 ) {
  FTFGamma0 = Gamma0;
}

// Set parameters of elastic scattering
inline void G4FTFParameters::SetAvaragePt2ofElasticScattering( const G4double aPt2 ) {
  AvaragePt2ofElasticScattering = aPt2;
}

// Set parameters of excitations

inline void G4FTFParameters::SetParams( const G4int ProcN,
                                        const G4double A1, const G4double B1, const G4double A2,
                                        const G4double B2, const G4double A3, const G4double Atop,
                                        const  G4double Ymin ) {
  ProcParams[ProcN][0] =   A1; ProcParams[ProcN][1] =  B1;
  ProcParams[ProcN][2] =   A2; ProcParams[ProcN][3] =  B2;
  ProcParams[ProcN][4] =   A3;
  ProcParams[ProcN][5] = Atop; ProcParams[ProcN][6] = Ymin;
}

inline void G4FTFParameters::SetDeltaProbAtQuarkExchange( const G4double aValue ) {
  DeltaProbAtQuarkExchange = aValue;
}

inline void G4FTFParameters::SetProbOfSameQuarkExchange( const G4double aValue ) {
  ProbOfSameQuarkExchange = aValue;
}

inline void G4FTFParameters::SetProjMinDiffMass( const G4double aValue ) {
  ProjMinDiffMass = aValue*CLHEP::GeV;
}

inline void G4FTFParameters::SetProjMinNonDiffMass( const G4double aValue ) {
  ProjMinNonDiffMass = aValue*CLHEP::GeV;
}

inline void G4FTFParameters::SetTarMinDiffMass( const G4double aValue ) {
  TarMinDiffMass = aValue*CLHEP::GeV;
}

inline void G4FTFParameters::SetTarMinNonDiffMass( const G4double aValue ) {
  TarMinNonDiffMass = aValue*CLHEP::GeV;
}

inline void G4FTFParameters::SetAveragePt2( const G4double aValue ) {
  AveragePt2 = aValue*CLHEP::GeV*CLHEP::GeV;
}

inline void G4FTFParameters::SetProbLogDistrPrD( const G4double aValue ) {
  ProbLogDistrPrD = aValue;
}

inline void G4FTFParameters::SetProbLogDistr( const G4double aValue ) {
  ProbLogDistr = aValue;
}

// Set parameters of a string kink

inline void G4FTFParameters::SetPt2Kink( const G4double aValue ) {
  Pt2kink = aValue;
}

inline void G4FTFParameters::SetQuarkProbabilitiesAtGluonSplitUp( const G4double Puubar, 
                                                                  const G4double Pddbar,
                                                                  const G4double Pssbar ) {
  QuarkProbabilitiesAtGluonSplitUp.push_back( Puubar ); 
  QuarkProbabilitiesAtGluonSplitUp.push_back( Puubar + Pddbar );
  QuarkProbabilitiesAtGluonSplitUp.push_back( Puubar + Pddbar + Pssbar );
}

// Set parameters of nuclear destruction
inline void G4FTFParameters::SetMaxNumberOfCollisions( const G4double Plab, 
                                                       const G4double Pbound ) {
  if ( Plab > Pbound ) {
    MaxNumberOfCollisions = Plab/Pbound;
    SetProbOfInteraction( -1.0 );
  } else {
    //MaxNumberOfCollisions = -1.0;
    //SetProbOfInteraction( G4Exp( 0.25*(Plab-Pbound) ) );
    MaxNumberOfCollisions = 1;
    SetProbOfInteraction( -1.0 );
  }
}

inline void G4FTFParameters::SetProbOfInteraction( const G4double aValue ) {
  ProbOfInelInteraction = aValue;
}

inline void G4FTFParameters::SetCofNuclearDestructionPr( const G4double aValue ) {
  CofNuclearDestructionPr = aValue;
}

inline void G4FTFParameters::SetCofNuclearDestruction( const G4double aValue ) {
  CofNuclearDestruction = aValue;
}

inline void G4FTFParameters::SetR2ofNuclearDestruction( const G4double aValue ) {
  R2ofNuclearDestruction = aValue;
}

inline void G4FTFParameters::SetExcitationEnergyPerWoundedNucleon( const G4double aValue ) {
  ExcitationEnergyPerWoundedNucleon = aValue;
}

inline void G4FTFParameters::SetDofNuclearDestruction( const G4double aValue ) {
  DofNuclearDestruction = aValue;
}

inline void G4FTFParameters::SetPt2ofNuclearDestruction( const G4double aValue ) {
  Pt2ofNuclearDestruction = aValue;
}

inline void G4FTFParameters::SetMaxPt2ofNuclearDestruction( const G4double aValue ) {
  MaxPt2ofNuclearDestruction = aValue;
}

// Get geometrical parameteres
inline G4double G4FTFParameters::GetTotalCrossSection() {
  return FTFXtotal;
}

inline G4double G4FTFParameters::GetElasticCrossSection() {
  return FTFXelastic;
}

inline G4double G4FTFParameters::GetInelasticCrossSection() {
  return FTFXinelastic;
}

inline G4double G4FTFParameters::GetSlope() {
  return FTFSlope;
}

inline G4double G4FTFParameters::GetProbabilityOfInteraction( const G4double impactsquare ) {
  if ( RadiusOfHNinteractions2 > impactsquare ) {
    return 1.0;
  } else {
    return 0.0;
  }
} 

inline G4double G4FTFParameters::GetProbabilityOfElasticScatt() {
  return ProbabilityOfElasticScatt;
}

inline G4double G4FTFParameters::GetInelasticProbability( const G4double impactsquare ) {
  G4double Gamma = GammaElastic( impactsquare );
  return 2*Gamma - Gamma*Gamma;
}

inline G4double G4FTFParameters::GetProbabilityOfAnnihilation() {
  return ProbabilityOfAnnihilation;
} 

// Get parameters of elastic scattering
inline G4double G4FTFParameters::GetAvaragePt2ofElasticScattering() {
  return AvaragePt2ofElasticScattering;
}

// Get parameters of excitations

inline G4double G4FTFParameters::GetDeltaProbAtQuarkExchange() { 
  return DeltaProbAtQuarkExchange;
}

inline G4double G4FTFParameters::GetProbOfSameQuarkExchange() {
  return ProbOfSameQuarkExchange;
}

inline G4double G4FTFParameters::GetProjMinDiffMass() {
  return ProjMinDiffMass;
}

inline G4double G4FTFParameters::GetProjMinNonDiffMass() {
  return ProjMinNonDiffMass;
}

inline G4double G4FTFParameters::GetTarMinDiffMass() {
  return TarMinDiffMass;
}

inline G4double G4FTFParameters::GetTarMinNonDiffMass() {
  return TarMinNonDiffMass;
}

inline G4double G4FTFParameters::GetAveragePt2() {
  return AveragePt2;
}

inline G4double G4FTFParameters::GetProbLogDistrPrD() {
  return ProbLogDistrPrD;
}

inline G4double G4FTFParameters::GetProbLogDistr() {
  return ProbLogDistr;
}

// Get parameters of a string kink

inline G4double G4FTFParameters::GetPt2Kink() {
  return Pt2kink;
}

inline std::vector< G4double > G4FTFParameters::GetQuarkProbabilitiesAtGluonSplitUp() {
  return QuarkProbabilitiesAtGluonSplitUp;
}

// Get parameters of nuclear destruction

inline G4double G4FTFParameters::GetMaxNumberOfCollisions() {
  return MaxNumberOfCollisions;
}

inline G4double G4FTFParameters::GetProbOfInteraction() {
  return ProbOfInelInteraction;
}

inline G4double G4FTFParameters::GetCofNuclearDestructionPr() {
  return CofNuclearDestructionPr;
}

inline G4double G4FTFParameters::GetCofNuclearDestruction() {
  return CofNuclearDestruction;
}

inline G4double G4FTFParameters::GetR2ofNuclearDestruction() {
  return R2ofNuclearDestruction;
}

inline G4double G4FTFParameters::GetExcitationEnergyPerWoundedNucleon() {
  return ExcitationEnergyPerWoundedNucleon;
}

inline G4double G4FTFParameters::GetDofNuclearDestruction() {
  return DofNuclearDestruction;
}

inline G4double G4FTFParameters::GetPt2ofNuclearDestruction() {
  return Pt2ofNuclearDestruction;
}

inline G4double G4FTFParameters::GetMaxPt2ofNuclearDestruction() {
  return MaxPt2ofNuclearDestruction;
}

#endif

