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
//---------------------------------------------------------------------------
//
// ClassName:      G4FTFTunings
//
// Author:         2022 Alberto Ribon
//
// Description:    Singleton to keep sets of parameters, called "tunes",
//                 for the FTF model.
//
//                 Please NOTE that, as of now (Fall 2022) ONLY ONE tune
//                 can be selected/applied; attempt to select multiple tunes
//                 will not results in any error messages, however further
//                 down the workflow only the FIRST of the activated tunes
//                 will be used.
//
//                 To use one of the tunes of this class, there is no need to
//                 change anything in this class, and use instead one of the
//                 following two UI commands, before initialization:
//                     /process/had/models/ftf/selectTuneByIndex integerIndex
//                 or  /process/had/models/ftf/selectTuneByName  stringName
//                 for instance:
//                     /process/had/models/ftf/selectTuneByIndex 1
//                     or
//                     /process/had/models/ftf/selectTuneByIndex 2
//                     or
//                     /process/had/models/ftf/selectTuneByIndex 3
//                 or
//                     /process/had/models/ftf/selectTuneByName baryon-tune2022-v0
//                     or
//                     /process/had/models/ftf/selectTuneByName pion-tune2022-v0
//                     or
//                     /process/had/models/ftf/selectTuneByName combined-tune2022-v0
//
//                 If you want to create a new tune, then you need to modify
//                 this class as follows: look for the first "dummy" tune
//                 available; if you find it, then specify its name in the
//                 std::array fNameOfTunes and the values of the parameters
//                 in the methods: G4FTFParamCollection::SetTuneN()
//                                 G4FTFParamCollBaryonProj::SetTuneN()
//                                 G4FTFParamCollMesonProj::SetTuneN()
//                                 G4FTFParamCollPionProj::SetTuneN
//                 Note that you need to set explicitly only the parameters
//                 with non-default values - all the others inherit the
//                 corresponding default values.
//                 If you don't find available "dummy" tune, then you need
//                 to increase by (at least) 1 the number of tunes, and add
//                 the corresponding "SetTuneN()" methods in the 4 classes
//                   G4FTFParamCollection, G4FTFParamCollBaryonProj,
//                   G4FTFParamCollMesonProj, G4FTFParamCollPionProj
//
//                 In order to explore some variations of FTF parameters
//                 (for instance to find out a new tune), please select
//                 (via UI command, as explained above) the existing tune
//                 from which you want to start with as "baseline", and
//                 then set the values of the parameters you want to change
//                 via the following C++ code (to used before initialization):
//                   G4HadronicDeveloperParameters::GetInstance()->Set(...)
//
//                 Note: in its current, first version, of this class,
//                       any FTF tune is applied "globally", i.e. for all
//                       projectile hadrons and regardless of their kinetic
//                       energy.
//                       In future versions, we might try to have tunes that
//                       are meant for specific projectile type and/or for
//                       intervals of kinetic energy (e.g. low-energy,
//                       medium-energy, high-energy).
//
//                 Note: a few classes (written by Julia Yarba) used only in
//                       G4FTFParameters, related to the set of parameters of
//                       the FTF models, have been moved from the header and
//                       source files of the class G4FTFParameters to this
//                       (G4FTFTunings) class, with minimal modifications.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4FTFTunings_h
#define G4FTFTunings_h 1

#include "globals.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include <array>

class G4ParticleDefinition;
class G4FTFTuningsMessenger;


class G4FTFTunings {
  public:
  
    static G4FTFTunings* Instance();
    ~G4FTFTunings();

    inline G4String GetTuneName( const G4int index ) const;
    // Returns the name of the specified tune (via its index).
    // Note that the name of the tune cannot be changed
    // (i.e. there is no corresponding "Set" method).
  
    inline G4int GetTuneApplicabilityState( const G4int index ) const;
    void SetTuneApplicabilityState( const G4int index, const G4int state );
    // Get/Set methods for the "applicability state" of the specified tune
    // (via its index). For the time being, there are only two states:
    // 0: switched off; 1: switched on.
  
    G4int GetIndexTune( const G4ParticleDefinition* particleDef, const G4double ekin ) const;
    // Based on the projectile type and its kinetic energy (from the input arguments),
    // this method returns the index of the tune which should be used.
    // For the time being, it returns the first alternative tune which is switched on,
    // else returns 0 which corresponds to the default set of parameters.
    // Note: this is the key method that needs to be revised if we decide to have
    //       different tunes according to projectile type and/or projectile energy range.
  
    static const G4int sNumberOfTunes = 10;
    // Number of tunes: must be >= 1, with the first one (i.e. with index = 0)
    // which corresponds to the default set of parameters.
    // For the time being, we set it to 10 : the second one (index = 1) is a
    // realistic alternative tune, whereas all the remaining 8 are "dummy" tunes,
    // i.e. the same as the default set of parameters. These are meant to be
    // replaced in the future with other, realistic alternative tunes.
    // Note: below, for the names and "applicability" status of tunes we use
    //       std::array - instead of std::vector - because the number of tunes
    //       do not change dynamically during a run, and, moreover, we expect
    //       quite a small number of them (just a few).
  
  private:
  
    G4FTFTunings();
    G4bool IsLocked() const;
  
    static G4FTFTunings* sInstance;

    G4FTFTuningsMessenger* fMessenger;

    const std::array< G4String, sNumberOfTunes > fNameOfTunes = { {
      "default",        // 0th tuning: default set
      "baryon-tune2022-v0",    // 1st tuning: Julia Yarba's presentation on 20-Jul-2022
      "pion-tune2022-v0",      // 2nd tuning: Julia Yarba's presentations on 26-Sept-2022 and 19-Oct-2022
      "combined-tune2022-v0",  // 3rd tuning: combo of the 1st and 2nd tuning 
      "fourth-dummy",   // 4th tuning: dummy
      "fifth-dummy",    // 5th tuning: dummy
      "sixth-dummy",    // 6th tuning: dummy
      "seventh-dummy",  // 7th tuning: dummy
      "eighth-dummy",   // 8th tuning: dummy
      "nineth-dummy"    // 9th tuning: dummy
    } };
    // The names of tunes can be useful for debugging.
  
    std::array< G4int, sNumberOfTunes > fApplicabilityOfTunes = { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    // Each tune has an integer that specifies its applicability.
    // For the time being, there only two values:
    //   0 : tune is switched off (i.e. not applicable);
    //   1 : tune is switched on (i.e. applicable).
    // Later on, it can be extended to indicate whether it is applicable to specific
    // projectile hadrons (e.g. protons, pions, etc.), and/or for specific energy ranges
    // (e.g. low-energy, medium-energy, high-energy - with energy thresholds to be
    // defined in this class).
    // The initial values can be changed (either via C++ interface or via UI command)
    // before initialization.
  
    //const G4double fLowEnergyThreshold  =  5.0*CLHEP::GeV;
    //const G4double fHighEnergyThreshold = 20.0*CLHEP::GeV;
    // These constants can be used, later on, to have different tunes
    // according to the energy of the projectile hadron (e.g. one set for
    // low energy, one set for middle energy, and one for high energy).
};


inline G4String G4FTFTunings::GetTuneName( const G4int index ) const {
  if ( index < 0 || index >= sNumberOfTunes ) return G4String();
  return fNameOfTunes[index];
}


inline G4int G4FTFTunings::GetTuneApplicabilityState( const G4int index ) const {
  if ( index < 0 || index >= sNumberOfTunes ) return 0;  // Switched off
  return fApplicabilityOfTunes[index];
}


//============================================================================

// Classes below have been created by Julia Yarba and were originally placed
// in the G4FTFParameters.{hh,cc} files ; some minimal changes and extensions
// have been included.


class G4FTFParamCollection {
  // NOTE: the settings are different for:
  //       * baryons projectile
  //       * anti-baryons projectile
  //       * pions (chg or pi0) projectile
  //       * kaons projectile (pdg = +/-321, 311, 130, or 310)
  //       * "undefined" projectile - nucleon assumed
  public:

    // Set-up the tune specified in the input argument, only if that tune is switched on.
    virtual void SetTune( const G4int tuneIndex );
  
    virtual void SetTune1();  // Set-up the 1st tune
    virtual void SetTune2();  // Set-up the 2nd tune
    virtual void SetTune3();  // Set-up the 3rd tune
    virtual void SetTune4();  // Set-up the 4th tune
    virtual void SetTune5();  // Set-up the 5th tune
    virtual void SetTune6();  // Set-up the 6th tune
    virtual void SetTune7();  // Set-up the 7th tune
    virtual void SetTune8();  // Set-up the 8th tune
    virtual void SetTune9();  // Set-up the 9th tune
    //...
  
    virtual ~G4FTFParamCollection() {}

    // parameters of excitation
    // Proc=0 --> Qexchg w/o excitation
    double GetProc0A1()   const  { return fProc0A1; }
    double GetProc0B1()   const  { return fProc0B1; }
    double GetProc0A2()   const  { return fProc0A2; }
    double GetProc0B2()   const  { return fProc0B2; }
    double GetProc0A3()   const  { return fProc0A3; }
    double GetProc0Atop() const  { return fProc0Atop; }
    double GetProc0Ymin() const  { return fProc0Ymin; }
    // Proc=1 --> Qexchg w/excitation
    double GetProc1A1()   const  { return fProc1A1; }
    double GetProc1B1()   const  { return fProc1B1; }
    double GetProc1A2()   const  { return fProc1A2; }
    double GetProc1B2()   const  { return fProc1B2; }
    double GetProc1A3()   const  { return fProc1A3; }
    double GetProc1Atop() const  { return fProc1Atop; }
    double GetProc1Ymin() const  { return fProc1Ymin; }
    // Proc=2 & Proc=3 in case ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )
    // Update: Proc=2 & Proc=3 in case ( AbsProjectileBaryonNumber > 10 ||  NumberOfTargetNucleons > 10 )
    // (diffraction dissociation)
    // Other parameters have a complex form for baryon projectile
    // although they're just numbers for e.g. pions projectile
    // Proc=2 --> Projectile diffraction
    double GetProc2A1()   const  { return fProc2A1; }
    double GetProc2B1()   const  { return fProc2B1; }
    double GetProc2A2()   const  { return fProc2A2; }
    double GetProc2B2()   const  { return fProc2B2; }
    double GetProc2A3()   const  { return fProc2A3; }
    double GetProc2Atop() const  { return fProc2Atop; }
    double GetProc2Ymin() const  { return fProc2Ymin; }
    // Proc=3 --> Target diffraction
    double GetProc3A1()   const  { return fProc3A1; }
    double GetProc3B1()   const  { return fProc3B1; }
    double GetProc3A2()   const  { return fProc3A2; }
    double GetProc3B2()   const  { return fProc3B2; }
    double GetProc3A3()   const  { return fProc3A3; }
    double GetProc3Atop() const  { return fProc3Atop; }
    double GetProc3Ymin() const  { return fProc3Ymin; }
    bool   IsProjDiffDissociation() const { return fProjDiffDissociation; }
    bool   IsTgtDiffDissociation()  const { return fTgtDiffDissociation; }
    // Proc=4 --> Qexchg "w/additional multiplier" in excitation
    double GetProc4A1()   const  { return fProc4A1; }
    double GetProc4B1()   const  { return fProc4B1; }
    double GetProc4A2()   const  { return fProc4A2; }
    double GetProc4B2()   const  { return fProc4B2; }
    double GetProc4A3()   const  { return fProc4A3; }
    double GetProc4Atop() const  { return fProc4Atop; }
    double GetProc4Ymin() const  { return fProc4Ymin; }
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
    // separately for baryons, mesons, etc.
    double GetR2ofNuclearDestruct()         const { return fR2ofNuclearDestruct; }
    double GetExciEnergyPerWoundedNucleon() const { return fExciEnergyPerWoundedNucleon; }
    double GetDofNuclearDestruct()          const { return fDofNuclearDestruct; } 
    double GetMaxPt2ofNuclearDestruct()     const { return fMaxPt2ofNuclearDestruct; }
   
  protected:

    G4FTFParamCollection();

    // parameters of excitation
    // these are for Inelastic interactions, i.e. Xinelastic=(Xtotal-Xelastix)>0.
    // for elastic, all the A's & B's, Atop & Ymin are zeros
    // general formula: Pp = A1*exp(B1*Y) + A2*exp(B2*Y) + A3
    // but if Y<Ymin, then Pp=max(0.,Atop)
    // for details, see also G4FTFParameters::GetProcProb( ProcN, y )
    // Proc=0 --> Qexchg w/o excitation
    double fProc0A1;
    double fProc0B1;
    double fProc0A2;
    double fProc0B2;
    double fProc0A3;
    double fProc0Atop;
    double fProc0Ymin;
    // Proc=1 --> Qexchg w/excitation
    double fProc1A1;
    double fProc1B1;
    double fProc1A2;
    double fProc1B2;
    double fProc1A3;
    double fProc1Atop;
    double fProc1Ymin;
    // NOTE: Proc #2 & 3 are projectile & target diffraction
    //       they have more complex definition of A1 & A2 
    //       for *baryons* although they're just numbers for pions
    //       (example for baryons below)
    // SetParams( 2, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Projectile diffraction
    // SetParams( 3, 6.0/Xinel, 0.0 ,-6.0/Xinel*16.28, 3.0 , 0.0, 0.0  ,     0.93);// Target diffraction
    //
    // Also, for ( AbsProjectileBaryonNumber > 1 ||  NumberOfTargetNucleons > 1 )
    // projectile and/or target diffraction (dissociation) may be switched ON/OFF 
    bool fProjDiffDissociation;
    bool fTgtDiffDissociation;
    // Proc=2 --> Projectile diffraction
    double fProc2A1; 
    double fProc2B1; 
    double fProc2A2; 
    double fProc2B2; 
    double fProc2A3; 
    double fProc2Atop; 
    double fProc2Ymin; 
    // Proc=3 --> Target diffraction
    double fProc3A1; 
    double fProc3B1; 
    double fProc3A2; 
    double fProc3B2; 
    double fProc3A3; 
    double fProc3Atop; 
    double fProc3Ymin; 
    // Proc=4 --> Qexchg w/additional multiplier in excitation  
    double fProc4A1;
    double fProc4B1;
    double fProc4A2;
    double fProc4B2;
    double fProc4A3;
    double fProc4Atop;
    double fProc4Ymin;
    // parameters of participating baryon excitation 
    // NOTE: baryon or HADRON ???
    // NOTE: this parameters (as C++ class data members) are used for all types of hadrons
    //       but the values for a specific group of particles can be are different from
    //       another group of particles
    //       the defaults listed under coments are for baryons, 
    //       and they may be different or the same for other hadrons (e.g. mesons)
    double fDeltaProbAtQuarkExchange;
    double fProbOfSameQuarkExchange;
    double fProjMinDiffMass;
    double fProjMinNonDiffMass;
    double fTgtMinDiffMass;
    double fTgtMinNonDiffMass;
    double fAveragePt2;
    double fProbLogDistrPrD;
    double fProbLogDistr;
    // parameters of nuclear distruction 
    // NOTE (JVY): there're 3 cases here:
    //             * baryon projectile
    //             * anti-baryon projectile
    //             * meson projectile
    // double fBaryonMaxNumberOfCollisions; // D=2.
    // void SetBaryonProbOfInteraction( const double ); // ??? this is prob. of inelastic interaction 
                                                        //     that is set internally based on certain conditions...
    // general (i.e. for used for baryons,anti-baryons, and mesons)
    // NOTE: these parameters have stayed THE SAME for quite a while 
    double fNuclearProjDestructP1;
    bool   fNuclearProjDestructP1_NBRNDEP;
    double fNuclearTgtDestructP1;
    bool   fNuclearTgtDestructP1_ADEP;
    double fNuclearProjDestructP2;
    double fNuclearProjDestructP3;
    double fNuclearTgtDestructP2;
    double fNuclearTgtDestructP3;
    //
    double fPt2NuclearDestructP1;
    double fPt2NuclearDestructP2;
    double fPt2NuclearDestructP3;
    double fPt2NuclearDestructP4;
    // baryons... well, in fact also mesons...
    double fR2ofNuclearDestruct;
    double fExciEnergyPerWoundedNucleon;
    double fDofNuclearDestruct;
    double fMaxPt2ofNuclearDestruct;
};


class G4FTFParamCollBaryonProj : public G4FTFParamCollection {
  public:
    G4FTFParamCollBaryonProj();
  
    virtual void SetTune1() override;  // Set-up the baryon part of the 1st tune
    virtual void SetTune2() override;  // Set-up the baryon part of the 2nd tune
    virtual void SetTune3() override;  // Set-up the baryon part of the 3rd tune
    virtual void SetTune4() override;  // Set-up the baryon part of the 4th tune
    virtual void SetTune5() override;  // Set-up the baryon part of the 5th tune
    virtual void SetTune6() override;  // Set-up the baryon part of the 6th tune
    virtual void SetTune7() override;  // Set-up the baryon part of the 7th tune
    virtual void SetTune8() override;  // Set-up the baryon part of the 8th tune
    virtual void SetTune9() override;  // Set-up the baryon part of the 9th tune
    //...
};


class G4FTFParamCollMesonProj : public G4FTFParamCollection {
  public:
    G4FTFParamCollMesonProj();

    virtual void SetTune1() override;  // Set-up the meson part of the 1st tune
    virtual void SetTune2() override;  // Set-up the meson part of the 2nd tune
    virtual void SetTune3() override;  // Set-up the meson part of the 3rd tune
    virtual void SetTune4() override;  // Set-up the meson part of the 4th tune
    virtual void SetTune5() override;  // Set-up the meson part of the 5th tune
    virtual void SetTune6() override;  // Set-up the meson part of the 6th tune
    virtual void SetTune7() override;  // Set-up the meson part of the 7th tune
    virtual void SetTune8() override;  // Set-up the meson part of the 8th tune
    virtual void SetTune9() override;  // Set-up the meson part of the 9th tune
    //...  
};


class G4FTFParamCollPionProj : public G4FTFParamCollMesonProj {
  public:    
    G4FTFParamCollPionProj();

    virtual void SetTune1() override;  // Set-up the pion part of the 1st tune
    virtual void SetTune2() override;  // Set-up the pion part of the 2nd tune
    virtual void SetTune3() override;  // Set-up the pion part of the 3rd tune
    virtual void SetTune4() override;  // Set-up the pion part of the 4th tune
    virtual void SetTune5() override;  // Set-up the pion part of the 5th tune
    virtual void SetTune6() override;  // Set-up the pion part of the 6th tune
    virtual void SetTune7() override;  // Set-up the pion part of the 7th tune
    virtual void SetTune8() override;  // Set-up the pion part of the 8th tune
    virtual void SetTune9() override;  // Set-up the pion part of the 9th tune
    //...
};

#endif
