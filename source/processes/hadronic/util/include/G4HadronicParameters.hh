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
// ClassName:      G4HadronicParameters
//
// Author:         2018 Alberto Ribon
//
// Description:    Singleton to keep global hadronic parameters.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronicParameters_h
#define G4HadronicParameters_h 1

#include "globals.hh"

class G4HadronicParametersMessenger;


class G4HadronicParameters {
  public:

    static G4HadronicParameters* Instance();
    ~G4HadronicParameters();

    inline G4double GetMaxEnergy() const;
    void SetMaxEnergy( const G4double val );
    // Getter/Setter for the upper limit for Geant4 hadronic physics, for any application.
    // Any hadronic model, physics list builder and constructor should use this method
    // instead of putting an arbitrary value in the code.
    // Any application which tries to use hadronic physics for an energy higher than this limit
    // will get a run-time crash, because no model is found.

    inline G4double GetMinEnergyTransitionFTF_Cascade() const;
    inline G4double GetMaxEnergyTransitionFTF_Cascade() const;
    void SetMinEnergyTransitionFTF_Cascade( const G4double val );
    void SetMaxEnergyTransitionFTF_Cascade( const G4double val );
    // Getter/Setter of the recommended energy limits, for physics lists, of the
    // transition region between the Fritiof (FTF) string model and the
    // intranuclear cascade model, either Bertini (BERT) or Binary (BIC). 

    inline G4double GetMinEnergyTransitionQGS_FTF() const;
    inline G4double GetMaxEnergyTransitionQGS_FTF() const;
    void SetMinEnergyTransitionQGS_FTF( const G4double val );
    void SetMaxEnergyTransitionQGS_FTF( const G4double val );
    // Getter/Setter of the recommended energy limits, for physics lists, of the
    // transition region between the two strings models - the Quark Gluon String (QGS)
    // model and the Fritiof (FTF) model.

    inline G4double GetMinEnergyINCLXX_Pbar() const;
    inline G4double GetMaxEnergyINCLXX_Pbar() const;
    void SetMinEnergyINCLXX_Pbar( const G4double val );
    void SetMaxEnergyINCLXX_Pbar( const G4double val );
    // Getter/Setter of the recommended energy limits, for physics lists, of the
    // intranuclear cascade model INCLXX, for pbar interaction. 

    inline G4double EnergyThresholdForHeavyHadrons() const;
    void SetEnergyThresholdForHeavyHadrons( G4double val );
    // If max kinetic energy is below this limit, then EM and hadronic physics are not 
    // instantiated for hyperons, anti-hyperons, anti light ions, b-, c- particles.

    inline G4double XSFactorNucleonInelastic() const;
    void SetXSFactorNucleonInelastic( G4double val );
    inline G4double XSFactorNucleonElastic() const;
    void SetXSFactorNucleonElastic( G4double val );
    // Cross section factor for protons and neutrons.

    inline G4double XSFactorPionInelastic() const;
    void SetXSFactorPionInelastic( G4double val );
    inline G4double XSFactorPionElastic() const;
    void SetXSFactorPionElastic( G4double val );
    // Cross section factor for pions.

    inline G4double XSFactorHadronInelastic() const;
    void SetXSFactorHadronInelastic( G4double val );
    inline G4double XSFactorHadronElastic() const;
    void SetXSFactorHadronElastic( G4double val );
    // Cross section factor for other hadrons and ions.

    inline G4double XSFactorEM() const;
    void SetXSFactorEM( G4double val );
    // Cross section factor for gamma and leptons.

    inline G4bool EnableBCParticles() const;
    void SetEnableBCParticles( G4bool val );
    // Baryons and mesons with c- and b- quarks may be enabled/disabled.
    // This flag is used both by EM and hadronic physics constructors.

    inline G4bool EnableHyperNuclei() const;
    void SetEnableHyperNuclei( G4bool val );
    // Light hyper-nuclei may be enabled/disabled.
    // This flag is used both by EM and hadronic physics constructors.

    inline G4bool ApplyFactorXS() const;
    void SetApplyFactorXS( G4bool val );
    // Flag enabling cross section factor definition.

    inline G4int GetVerboseLevel() const;
    void SetVerboseLevel( const G4int val );
    // Getter/Setter of the general verbosity level for hadronics.
  
    inline G4bool EnableCRCoalescence() const;
    void SetEnableCRCoalescence( G4bool val );
    // Boolean switch that allows to apply the Cosmic Ray (CR) coalescence algorithm
    // to the secondaries produced by a string model. By default it is disabled.

    inline G4bool EnableIntegralInelasticXS() const;
    inline G4bool EnableIntegralElasticXS() const;
    void SetEnableIntegralInelasticXS( G4bool val );
    void SetEnableIntegralElasticXS( G4bool val );
    // Enable/disable integral method for main types of hadrons.
  
    inline G4bool EnableDiffDissociationForBGreater10() const;
    // For nucleon-hadron interactions, it's not decided what to do with diffraction
    // dissociation. For the moment, they are turned off. This option allows it to
    // be turned back on. Applies to Baryon Number > 10 or # target nucleons > 10.
    void SetEnableDiffDissociationForBGreater10(G4bool val);

    inline G4bool EnableCoherentChargeExchange() const;
    void SetEnableCoherentChargeExchange( G4bool val );
    // Coherent Charge exchange process may be enabled/disabled.

    inline G4bool EnableNeutronGeneralProcess() const;
    void SetEnableNeutronGeneralProcess( G4bool val );
    // Neutron general process may be enabled/disabled.

    inline G4double GetEPRelativeLevel() const;
    inline G4double GetEPAbsoluteLevel() const;
    inline G4int GetEPReportLevel() const;
    inline G4bool GetBinaryDebug() const;
    inline const G4String& GetDirPARTICLEXS() const;
    inline const G4String& GetPhysListDocDir() const;
    inline const G4String& GetPhysListName() const;
    // Access to environment variables.

    inline G4double GetNeutronKineticEnergyThresholdForSVT() const;
    void SetNeutronKineticEnergyThresholdForSVT( const G4double val );
    // Getter/Setter for the neutron kinetic energy threshold for 
    // applying the SVT (Sampling of the Velocity of the Target) algorithm.

    inline G4double GetTimeThresholdForRadioactiveDecay() const;
    void SetTimeThresholdForRadioactiveDecay( const G4double val );
    // Getter/Setter for the time threshold of radioactive decays
    // (i.e. radioactive decays that happen later than this value are ignored).

  private:

    G4HadronicParameters();

    G4bool IsLocked() const;

    static G4HadronicParameters* sInstance;

    G4HadronicParametersMessenger* fMessenger;

    G4double fMaxEnergy;
    G4double fMinEnergyTransitionFTF_Cascade;
    G4double fMaxEnergyTransitionFTF_Cascade;
    G4double fMinEnergyTransitionQGS_FTF;
    G4double fMaxEnergyTransitionQGS_FTF;
    G4double fMinEnergyINCLXX_Pbar;
    G4double fMaxEnergyINCLXX_Pbar;
    G4double fEnergyThresholdForHeavyHadrons;
    G4double fXSFactorNucleonInelastic = 1.0;
    G4double fXSFactorPionInelastic = 1.0;
    G4double fXSFactorHadronInelastic = 1.0;
    G4double fXSFactorNucleonElastic = 1.0;
    G4double fXSFactorPionElastic = 1.0;
    G4double fXSFactorHadronElastic = 1.0;
    G4double fXSFactorEM = 1.0;
    G4double fXSFactorLimit = 0.2;
    G4double fRelativeDiff = DBL_MAX;
    G4double fAbsoluteDiff = DBL_MAX;
    G4double fNeutronEkinThresholdForSVT = -1.0;
    G4double fTimeThresholdForRadioactiveDecays = -1.0;
    
    G4int fVerboseLevel = 1;
    G4int fReportLevel = 0;

    G4bool fEnableBC = false;
    G4bool fEnableHyperNuclei = false;
    G4bool fApplyFactorXS = false;
    G4bool fEnableCRCoalescence = false;
    G4bool fEnableIntegralInelasticXS = true;
    G4bool fEnableIntegralElasticXS = true;
    G4bool fEnableDiffDissociationForBGreater10 = false;
    G4bool fNeutronGeneral = false;
    G4bool fChargeExchange = false;
    G4bool fBinaryDebug = false;

    G4String fDirPARTICLEXS = "";
    G4String fPhysListDocDir = "";
    G4String fPhysListName = "";
};

inline G4double G4HadronicParameters::GetMaxEnergy() const { 
  return fMaxEnergy;
}

inline G4double G4HadronicParameters::GetMinEnergyTransitionFTF_Cascade() const { 
  return fMinEnergyTransitionFTF_Cascade;
}
inline G4double G4HadronicParameters::GetMaxEnergyTransitionFTF_Cascade() const { 
  return fMaxEnergyTransitionFTF_Cascade;
}

inline G4double G4HadronicParameters::GetMinEnergyTransitionQGS_FTF() const { 
  return fMinEnergyTransitionQGS_FTF;
}

inline G4double G4HadronicParameters::GetMaxEnergyTransitionQGS_FTF() const { 
  return fMaxEnergyTransitionQGS_FTF;
}

inline G4double G4HadronicParameters::GetMinEnergyINCLXX_Pbar() const { 
  return fMinEnergyINCLXX_Pbar;
}
inline G4double G4HadronicParameters::GetMaxEnergyINCLXX_Pbar() const { 
  return fMaxEnergyINCLXX_Pbar;
} 
  

inline G4double G4HadronicParameters::EnergyThresholdForHeavyHadrons() const {
  return fEnergyThresholdForHeavyHadrons;
}

inline G4double G4HadronicParameters::XSFactorNucleonInelastic() const {
  return fXSFactorNucleonInelastic;
}

inline G4double G4HadronicParameters::XSFactorNucleonElastic() const {
  return fXSFactorNucleonElastic;
}

inline G4double G4HadronicParameters::XSFactorPionInelastic() const {
  return fXSFactorPionInelastic;
}

inline G4double G4HadronicParameters::XSFactorPionElastic() const {
  return fXSFactorPionElastic;
}

inline G4double G4HadronicParameters::XSFactorHadronInelastic() const {
  return fXSFactorHadronInelastic;
}

inline G4double G4HadronicParameters::XSFactorHadronElastic() const {
  return fXSFactorHadronElastic;
}

inline G4double G4HadronicParameters::XSFactorEM() const {
  return fXSFactorEM;
}

inline G4int G4HadronicParameters::GetVerboseLevel() const { 
  return fVerboseLevel;
}

inline G4bool G4HadronicParameters::EnableBCParticles() const {
  return fEnableBC;
}

inline G4bool G4HadronicParameters::EnableHyperNuclei() const {
  return fEnableHyperNuclei;
}

inline G4bool G4HadronicParameters::ApplyFactorXS() const {
  return fApplyFactorXS;
}

inline G4bool G4HadronicParameters::EnableCRCoalescence() const {
  return fEnableCRCoalescence;
}

inline G4bool G4HadronicParameters::EnableIntegralInelasticXS() const {
  return fEnableIntegralInelasticXS;
}

inline G4bool G4HadronicParameters::EnableIntegralElasticXS() const {
  return fEnableIntegralElasticXS;
}

inline G4bool G4HadronicParameters::EnableDiffDissociationForBGreater10() const {
  return fEnableDiffDissociationForBGreater10;
}

inline G4bool G4HadronicParameters::EnableNeutronGeneralProcess() const {
  return fNeutronGeneral;
}

inline G4bool G4HadronicParameters::EnableCoherentChargeExchange() const {
  return fChargeExchange;
}

inline G4bool G4HadronicParameters::GetBinaryDebug() const {
  return fBinaryDebug;
}

inline G4double G4HadronicParameters::GetEPRelativeLevel() const {
  return fRelativeDiff;
}

inline G4double G4HadronicParameters::GetEPAbsoluteLevel() const {
  return fAbsoluteDiff;
}

inline G4int G4HadronicParameters::GetEPReportLevel() const {
  return fReportLevel;
}

inline const G4String& G4HadronicParameters::GetDirPARTICLEXS() const {
  return fDirPARTICLEXS;
}

inline const G4String& G4HadronicParameters::GetPhysListDocDir() const
{
  return fPhysListDocDir;
}

inline const G4String& G4HadronicParameters::GetPhysListName() const
{
  return fPhysListName;
}

inline G4double G4HadronicParameters::GetNeutronKineticEnergyThresholdForSVT() const { 
  return fNeutronEkinThresholdForSVT;
}

inline G4double G4HadronicParameters::GetTimeThresholdForRadioactiveDecay() const { 
  return fTimeThresholdForRadioactiveDecays;
}

#endif
