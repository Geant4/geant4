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
// $Id$
//
// GEANT4 physics class: G4QHadronInelasticDataSet -- header file
// Created by M. Kosov (Mikhail.Kossov@cern.ch) 11.11.09
//
// ------------------------------------------------------------------------
// Short description: G4hadr wrapper for CHIPS inelastic hA cross-sections.
// ------------------------------------------------------------------------

#ifndef G4QHadronInelasticDataSet_h
#define G4QHadronInelasticDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4SigmaZero.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiProton.hh"
#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"
//#include "G4Gamma.hh"
//#include "G4MuonPlus.hh"
//#include "G4MuonMinus.hh"
//#include "G4Electron.hh"
//#include "G4Positron.hh"
//#include "G4TauPlus.hh"
//#include "G4TauMinus.hh"
//#include "G4AntiNeutrinoE.hh"
//#include "G4NeutrinoE.hh"
//#include "G4AntiNeutrinoMu.hh"
//#include "G4NeutrinoMu.hh"
//#include "G4AntiNeutrinoTau.hh"
//#include "G4NeutrinoTau.hh"

#include "G4QIsotope.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QNeutronNuclearCrossSection.hh"  
#include "G4QPionMinusNuclearCrossSection.hh"
#include "G4QPionPlusNuclearCrossSection.hh"
#include "G4QKaonPlusNuclearCrossSection.hh"
#include "G4QKaonMinusNuclearCrossSection.hh"
#include "G4QKaonZeroNuclearCrossSection.hh"
#include "G4QHyperonNuclearCrossSection.hh"
#include "G4QHyperonPlusNuclearCrossSection.hh"
#include "G4QAntiBaryonPlusNuclearCrossSection.hh"
#include "G4QAntiBaryonNuclearCrossSection.hh"
//#include "G4QPhotonNuclearCrossSection.hh"
//#include "G4QElectronNuclearCrossSection.hh"
//#include "G4QMuonNuclearCrossSection.hh"
//#include "G4QTauNuclearCrossSection.hh"
//#include "G4QNuMuNuclearCrossSection.hh"
//#include "G4QANuMuNuclearCrossSection.hh"
//#include "G4QNuENuclearCrossSection.hh"
//#include "G4QANuENuclearCrossSection.hh"
//#include "G4QNuNuNuclearCrossSection.hh"
//#include "G4QANuANuNuclearCrossSection.hh"   


class G4QHadronInelasticDataSet : public G4VCrossSectionDataSet
{
  public:

    G4QHadronInelasticDataSet(const G4String& input_name = "CHIPSInelasticXS");
    ~G4QHadronInelasticDataSet() {};
    G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,    
                           const G4Element* elm = 0, const G4Material* mat = 0);
    G4double GetIsoCrossSection(const G4DynamicParticle* P, G4int Z, G4int A,
			              const G4Isotope*, const G4Element*, const G4Material*);
    void BuildPhysicsTable(const G4ParticleDefinition&){};
    void DumpPhysicsTable(const G4ParticleDefinition&) {};

    void Description() const;

  //private:

  //G4QIsotope* Isotopes;      // Pointer to the G4QIsotopes singleton
  //static std::vector <G4int> ElementZ;
  //// Z of the element(i) in theLastCalc
  //
  //static std::vector <std::vector<G4int>*> ElIsoN;
  //// N of isotope(j) of Element(i)
  //
  //static std::vector <std::vector<G4double>*> IsoProbInEl;
  //// SumProbabIsotopes in Element i
};

#endif
