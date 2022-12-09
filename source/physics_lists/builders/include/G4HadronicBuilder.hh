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
// Geant4 header G4HadronicBuilder
//
// Author V.Ivanchenko 14.05.2020
//
// Build hadronic physics 
//

#ifndef G4HadronicBuilder_h
#define G4HadronicBuilder_h 1

#include "globals.hh"
#include <vector>

class G4HadronicBuilder
{
private:

  // generic methods, Glauber-Gribov cross sections are used.
  // if the boolean "bert" is true (false) then Bertini cascade is (not) built;
  // when "bert" is false, then FTFP is used down to zero kinetic energy.
  
  static void BuildFTFP_BERT(const std::vector<G4int>& particleList, 
                             G4bool bert, const G4String& xsName);

  static void BuildFTFQGSP_BERT(const std::vector<G4int>& particleList, 
                                G4bool bert, const G4String& xsName);

  static void BuildQGSP_FTFP_BERT(const std::vector<G4int>& particleList, 
                                  G4bool bert, G4bool quasiElastic,
                                  const G4String& xsName);

public:

  // methods to build elastic and inelastic physics per particle category
  static void BuildElastic(const std::vector<G4int>& particleList);

  static void BuildHyperonsFTFP_BERT();

  static void BuildHyperonsFTFQGSP_BERT();

  static void BuildHyperonsQGSP_FTFP_BERT(G4bool quasiElastic);

  static void BuildKaonsFTFP_BERT();

  static void BuildKaonsFTFQGSP_BERT();

  static void BuildKaonsQGSP_FTFP_BERT(G4bool quasiElastic);

  static void BuildAntiLightIonsFTFP();

  //static void BuildAntiLightIonsQGSP_FTFP(G4bool quasiElastic);

  static void BuildBCHadronsFTFP_BERT();

  static void BuildBCHadronsFTFQGSP_BERT();

  static void BuildBCHadronsQGSP_FTFP_BERT(G4bool quasiElastic);

  // method to create some decays for heavy hadrons
  static void BuildDecayTableForBCHadrons();

  // methods for nuclear interactions of light hypernuclei and anti-hypernuclei projectiles
  // (note that: QGS cannot handle nuclei projectiles of any kind; INCLXX is currently
  // the only intra-nuclear cascade model capable of handling light hypernuclei - but not
  // light anti-hypernuclei - so only the two reference physics lists FTFP_BERT and 
  // FTFP_INCLXX are able to simulate nuclear interactions of light hypernuclei and
  //  anti-hypernuclei).
  static void BuildHyperNucleiFTFP_BERT();
  static void BuildHyperAntiNucleiFTFP_BERT();
  static void BuildHyperNucleiFTFP_INCLXX();
  static void BuildFTFP_INCLXX( const std::vector< G4int >& partList, const G4String& xsName );
};

#endif
