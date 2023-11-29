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
// ClassName:   G4HadronInelasticQBBC_ABLA
//
// Author: Alberto Ribon (CERN), April 2023
//
// Similar to the physics list constructor G4HadronInelasticQBBC_ABLA,
// except for the final-state of inelastic interactions of charged pions and
// nucleons in which ABLA nuclear de-excitation is utilized (instead of the
// usual Precompound/de-excitation).
// This is meant for testing purposes of the coupling between the hadronic
// string models (FTF and QGS) and ABLA (via G4GeneratorPrecompoundInterface),
// as well as of the coupling between intra-nuclear cascade models (BERT and BIC)
// and ABLA.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronInelasticQBBC_ABLA_h
#define G4HadronInelasticQBBC_ABLA_h 1

#include "globals.hh"
#include "G4VHadronPhysics.hh"


class G4HadronInelasticQBBC_ABLA : public G4VHadronPhysics {
  public: 
    G4HadronInelasticQBBC_ABLA( G4int ver = 1 );
  virtual ~G4HadronInelasticQBBC_ABLA() = default;
  void ConstructProcess() override;
  G4HadronInelasticQBBC_ABLA( G4HadronInelasticQBBC_ABLA& ) = delete;
  G4HadronInelasticQBBC_ABLA& operator=( const G4HadronInelasticQBBC_ABLA& right ) = delete;
};

#endif
