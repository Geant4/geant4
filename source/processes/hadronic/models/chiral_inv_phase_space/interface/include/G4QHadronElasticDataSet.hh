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
// $Id: G4QHadronElasticDataSet.hh,v 1.2 2010-02-16 07:53:05 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// GEANT4 physics class: G4QHadronElasticDataSet -- header file
// Created by M. Kosov (Mikhail.Kossov@cern.ch) 21.01.10
//
// ----------------------------------------------------------------------
// Short description: G4hadr wrapper for CHIPS elastic hA cross-sections.
// ----------------------------------------------------------------------

#ifndef G4QHadronElasticDataSet_h
#define G4QHadronElasticDataSet_h 1

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

#include "G4QIsotope.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"  
#include "G4QPionMinusElasticCrossSection.hh"
#include "G4QPionPlusElasticCrossSection.hh"
#include "G4QKaonPlusElasticCrossSection.hh"
#include "G4QKaonMinusElasticCrossSection.hh"
#include "G4QHyperonElasticCrossSection.hh"
#include "G4QHyperonPlusElasticCrossSection.hh"
#include "G4QAntiBaryonElasticCrossSection.hh"


class G4QHadronElasticDataSet : public G4VCrossSectionDataSet
{
public:

  G4QHadronElasticDataSet();
  ~G4QHadronElasticDataSet() {}
  G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Element* anElement);
  G4bool IsZAApplicable(const G4DynamicParticle* aParticle, G4double Z, G4double A);
  G4double GetCrossSection(const G4DynamicParticle* Part, const G4Element* El, G4double T);
  G4double GetIsoZACrossSection(const G4DynamicParticle* P,G4double Z,G4double A,G4double);
  void BuildPhysicsTable(const G4ParticleDefinition&){}
  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:

  G4QIsotope* Isotopes;                             // Pointer to the G4QIsotopes singleton
  static std::vector <G4int> ElementZ;                // Z of the element(i) in theLastCalc
  static std::vector <std::vector<G4int>*> ElIsoN;    // N of isotope(j) of Element(i)
  static std::vector <std::vector<G4double>*> IsoProbInEl;// SumProbabIsotopes in Element i
};

#endif
