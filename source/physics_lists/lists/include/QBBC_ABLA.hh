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
// ClassName:  QBBC_ABLA
//
// Author: Alberto Ribon (CERN), April 2023
//
// The new, experimental physics list QBBC_ABLA is similar to the reference
// physics list QBBC, except that for hadron inelastic the physics constructor
// G4HadronInelasticQBBC_ABLA is used (instead of G4HadronInelasticQBBC):
// in practice, QBBC_ABLA behaves as QBBC, with the only difference that for
// the final-state of nuclear inelastic interactions of charged pions and
// nucleons projectiles, the ABLA model (instead of the usual
// Precompound/de-excitation) is utilized for nuclear de-excitation.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef QBBC_ABLA_h
#define QBBC_ABLA_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"


class QBBC_ABLA : public G4VModularPhysicsList {
  public:
    explicit QBBC_ABLA( G4int ver = 1, const G4String& type = "QBBC_ABLA" );
    virtual ~QBBC_ABLA() = default;
    QBBC_ABLA( const QBBC_ABLA& ) = delete;
    QBBC_ABLA& operator=( const QBBC_ABLA& right ) = delete;
};

#endif
