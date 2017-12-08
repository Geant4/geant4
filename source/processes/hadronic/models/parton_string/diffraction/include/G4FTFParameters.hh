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
// $Id: G4FTFParameters.hh 107317 2017-11-08 16:25:57Z gcosmo $
// GEANT4 tag $Name:  $
//
#ifndef G4FTFParameters_h
#define G4FTFParameters_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ChipsComponentXS.hh"

#include "G4Exp.hh"

     // NOTE: the settings are different for:
     //       * baryons projectile
     //       * anti-baryons projectile
     //       * pions (chg or pi0) projectile
     //       * kaons projectile (pdg = +/-321, 311, 130, or 310)
     //       * "undefined" projectile - nucleon assumed

class G4FTFParamCollection {
  public:
    //dtor
    virtual ~G4FTFParamCollection() {}

    // parameters of excitation
    //
    // Proc=0 --> Qexchg w/o excitation
    //
    double GetProc0A1()   const  { return fProc0A1; }
    double GetProc0B1()   const  { return fProc0B1; }
    double GetProc0A2()   const  { return fProc0A2; }
    double GetProc0B2()   const  { return fProc0B2; }
    double GetProc0A3()   const  { return fProc0A3; }
    double GetProc0Atop() const  { return fProc0Atop; }
    double GetProc0Ymin() const  { return fProc0Ymin; }
    //
    // Proc=1 --> Qexchg w/excitation
    //
    // Proc=2 & Proc=3 for the case ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )
    // (diffraction dissociation)
    //
    bool   IsProjDiffDissociation() const { return fProjDiffDissociation; }
    bool   IsTgtDiffDissociation()  const { return fTgtDiffDissociation; }
    //
    double GetProc1A1()   const  { return fProc1A1; }
    double GetProc1B1()   const  { return fProc1B1; }
    double GetProc1A2()   const  { return fProc1A2; }
    double GetProc1B2()   const  { return fProc1B2; }
    double GetProc1A3()   const  { return fProc1A3; }
    double GetProc1Atop() const  { return fProc1Atop; }
    double GetProc1Ymin() const  { return fProc0Ymin; }
    //
    // Proc=4 --> Qexchg "w/additional multiplier" in excitation
    //
    double GetProc4A1()   const  { return fProc4A1; }
    double GetProc4B1()   const  { return fProc4B1; }
    double GetProc4A2()   const  { return fProc4A2; }
    double GetProc4B2()   const  { return fProc4B2; }
    double GetProc4A3()   const  { return fProc4A3; }
    double GetProc4Atop() const  { return fProc4Atop; }
    double GetProc4Ymin() const  { return fProc4Ymin; }
    //
    // 
    double GetDeltaProbAtQuarkExchange() const  { return fDeltaProbAtQuarkExchange; }  
    double GetProbOfSameQuarkExchange()  const  { return fProbOfSameQuarkExchange; }
    double GetProjMinDiffMass()          const  { return fProjMinDiffMass; }
    double GetProjMinNonDiffMass()       const  { return fProjMinNonDiffMass; }
    double GetTgtMinDiffMass()           const  { return fTgtMinDiffMass; }
    double GetTgtMinNonDiffMass()        const  { return fTgtMinNonDiffMass; }
    double GetAveragePt2()               const  { return fAveragePt2; }
    double GetProbLogDistrPrD()          const  { return fProbLogDistrPrD; }
    double GetProbLogDistr()             const  { return fProbLogDistr; }

    // NOTE (JVY): There is also the Pt2Kind parameter but for now it's set to 0., so we'll leave it aside
    // --> FIXME !!! --> void Get/SetBaryonMaxNumberOfCollisions( const double, const double ); // 1st is Plab, 2nd - D=2.
    // NOTE (JVY): These parameters are COMMON among various projectiles !!!
    //
    double GetNuclearProjDestructP1()    const { return fNuclearProjDestructP1; }
    bool   IsNuclearProjDestructP1_NBRNDEP() const { return fNuclearProjDestructP1_NBRNDEP; }
    double GetNuclearTgtDestructP1()     const { return fNuclearTgtDestructP1; }
    bool   IsNuclearTgtDestructP1_ADEP() const { return fNuclearTgtDestructP1_ADEP; }
    double GetNuclearProjDestructP2()    const { return fNuclearProjDestructP2; }
    double GetNuclearProjDestructP3()    const { return fNuclearProjDestructP3; }
    double GetNuclearTgtDestructP2()     const { return fNuclearTgtDestructP2; }
    double GetNuclearTgtDestructP3()     const { return fNuclearTgtDestructP3; }
    double GetPt2NuclearDestructP1()     const { return fPt2NuclearDestructP1; }
    double GetPt2NuclearDestructP2()     const { return fPt2NuclearDestructP2; }
    double GetPt2NuclearDestructP3()     const { return fPt2NuclearDestructP3; }
    double GetPt2NuclearDestructP4()     const { return fPt2NuclearDestructP4; }      
    //
    // separately for baryons, mesons, etc.
    //
    double GetR2ofNuclearDestruct()         const { return fR2ofNuclearDestruct; }
    double GetExciEnergyPerWoundedNucleon() const { return fExciEnergyPerWoundedNucleon; }
    double GetDofNuclearDestruct()          const { return fDofNuclearDestruct; } 
    double GetMaxPt2ofNuclearDestruct()     const { return fMaxPt2ofNuclearDestruct; }
   
  protected:
    // ctor
    G4FTFParamCollection();

    // parameters of excitation
    //
    //
    // these are for Inelastic interactions, i.e. Xinelastic=(Xtotal-Xelastix)>0.
    // for elastic, all the A's & B's, Atop & Ymin are zeros
    // general formula: Pp = A1*exp(B1*Y) + A2*exp(B2*Y) + A3
    // but if Y<Ymin, then Pp=max(0.,Atop)
    // for details, see also G4FTFParameters::GetProcProb( ProcN, y )
    //
    // Proc=0 --> Qexchg w/o excitation
    double fProc0A1; // D=13.71 
    double fProc0B1; // D=1.75
    double fProc0A2; // D=-30.69 (or -214.5 as in Doc ?)
    double fProc0B2; // D=3.     ( or 4. as in Doc ?)
    double fProc0A3; // D=0.
    double fProc0Atop; // D=1.   ( or 0.5 as in Doc ?)
    double fProc0Ymin; // D=0.93 (or 1.1 as in Doc ?)
    // Proc=1 --> Qexchg w/excitation
    double fProc1A1; // D=25.
    double fProc1B1; // D=1.
    double fProc1A2; // D=-50.34
    double fProc1B2; // D=1.5
    double fProc1A3; // D=0.
    double fProc1Atop; // D=0.
    double fProc1Ymin; // D=1.4
    //
    // NOTE: Proc #2 & 3 are projectile & target diffraction
    //       they have more complex definition of A1 & A2 
    //      (see around line 540 or so)
    // SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Projectile diffraction
    // SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Target diffraction
    //
    // Also, for ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )
    // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
    bool   fProjDiffDissociation;
    bool   fTgtDiffDissociation;
    //
    // Proc=4 --> Qexchg w/additional multiplier in excitation  
    double fProc4A1; // D=0.6 (or 1. as in Doc ?)
    double fProc4B1; // D=0.
    double fProc4A2; // D=-1.2 (or -2.01 as in Doc ?)
    double fProc4B2; // D=0.5 
    double fProc4A3; // D=0.
    double fProc4Atop; // D=0.
    double fProc4Ymin; // D=1.4
    //
    // parameters of participating baryon excitation
    //
    double fDeltaProbAtQuarkExchange; // D=0. 
    double fProbOfSameQuarkExchange;  // D=0. if A<=26, otherwise D=1.
    double fProjMinDiffMass;          // projectile, D=1.16GeV
    double fProjMinNonDiffMass;       // projectile, D=1.16GeV
    double fTgtMinDiffMass;           // target, D=1.16GeV
    double fTgtMinNonDiffMass;        // target, D=1.16GeV
    double fAveragePt2;               // D=0.3GeV**2 ( or 0.15 as in the Doc ???)
    double fProbLogDistrPrD;          // D=0.6 (or 0.3 ???)
    double fProbLogDistr;             // D=0.6 (or 0.3 ???)

    // parameters of nuclear distruction 
    //
    // NOTE (JVY): there're 3 cases here:
    //             * baryon projectile
    //             * anti-baryon projectile
    //             * meson projectile
    //
    // double fBaryonMaxNumberOfCollisions; // D=2.
    // void SetBaryonProbOfInteraction( const double ); // ??? this is prob. of inelastic interaction 
                                                        //     that is set internally based on certain conditions...
    // general (i.e. for used for baryons,anti-baryons, and mesons)
    // NOTE: these parameters have stayed THE SAME for quite a while 
    double fNuclearProjDestructP1; // D=0.00481 in 10.3.ref04 !!!
                                   // BUT !!! In 10.3.ref04 as well as in 10.2-seriesit's multiplied of AbsProjectileBaryonNumber
		                   // which somehow is 0 for the proton projectile (see in 10.3.ref04 around lines 130-140 In G4FTFParameters.cc).
				   // For the target destr. it's multipled by the number of target nucleons (12 for Carbon).
				   // In 10.3.p01 it's set to 1. FLAT OUT for both projectile & target, no multiplications, etc.
				   // Now, make default at 1.
    bool   fNuclearProjDestructP1_NBRNDEP;
    double fNuclearTgtDestructP1;  // Make D=1. as in 10.3.p01
    bool   fNuclearTgtDestructP1_ADEP;
    double fNuclearProjDestructP2; // D=4.0
    double fNuclearProjDestructP3; // D=2.1
    double fNuclearTgtDestructP2; // D=4.0
    double fNuclearTgtDestructP3; // D=2.1
    //
    double fPt2NuclearDestructP1; // D=0.035
    double fPt2NuclearDestructP2; // D=0.04
    double fPt2NuclearDestructP3; // D=4.0
    double fPt2NuclearDestructP4; // D=2.5 
    // baryons
    double fR2ofNuclearDestruct;         // D=1.5*fermi*fermi
    double fExciEnergyPerWoundedNucleon; // D=40MeV
    double fDofNuclearDestruct;          // D=0.3
    // NOTE: this parameter has changed from 1. to 9. between 10.2 and 10.4.ref04 !!!
    double fMaxPt2ofNuclearDestruct;     // D=9GeV**2

  private:
    void Reset();
};


class G4FTFParamCollBaryonProj : public G4FTFParamCollection {
  public:
    // ctor 
    G4FTFParamCollBaryonProj();
};


class G4FTFParameters {
  public:
    // G4FTFParameters( const G4ParticleDefinition* , G4int theA, G4int theZ, G4double s );
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
    //private: 
    // --->    G4FTFParameters();

    // Initial energy of hN interactions
    G4double FTFhNcmsEnergy;  // Initial hN CMS energy

    // hN cross section manager
    G4ChipsComponentXS* FTFxsManager;

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

    G4double DofNuclearDestruction;       // D for momentum sampling
    G4double Pt2ofNuclearDestruction;     // Pt2
    G4double MaxPt2ofNuclearDestruction;  // Max Pt2

  private:
  
    void Reset(); 

    // JVY, Oct. 31, 2017: encapsulates (current set of) parameters for the baryon projectile
    //
    G4FTFParamCollBaryonProj fParCollBaryonProj;

    // G4-MT changes
  private:
    static G4ThreadLocal bool chipsComponentXSisInitialized;
    static G4ThreadLocal G4ChipsComponentXS* chipsComponentXSinstance;
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

