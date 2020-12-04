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
// Author: Alberto Ribon
// Date:   October 2017
//
// Hadron physics for the new, experimental physics list FTFQGSP_BERT,
// with QGS fragmentation of strings, instead of the Lund string
// fragmentation. Note that the string excitation is still done with FTF,
// exactly as for FTFP_BERT.
// Given that it is an experimental, and perhaps temporary, new type of
// hadron physics, corresponding builders are not created and everything
// is implemented directly in this class.
//
// Modified:
// 19.10.2020 V.Ivanchenko use inheritance from G4HadronPhysicsFTFP_BERT
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTFQGSP_BERT_h
#define G4HadronPhysicsFTFQGSP_BERT_h 1

#include "G4HadronPhysicsFTFP_BERT.hh"

class G4HadronPhysicsFTFQGSP_BERT : public G4HadronPhysicsFTFP_BERT
{
  public: 
    G4HadronPhysicsFTFQGSP_BERT(G4int verbose =1);
    G4HadronPhysicsFTFQGSP_BERT(const G4String& name, G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsFTFQGSP_BERT();

    void ConstructProcess() override;

    // copy constructor and hide assignment operator
    G4HadronPhysicsFTFQGSP_BERT(G4HadronPhysicsFTFQGSP_BERT &) = delete;
    G4HadronPhysicsFTFQGSP_BERT & operator =
    (const G4HadronPhysicsFTFQGSP_BERT &right) = delete;

  protected:
    void DumpBanner() override;

};

#endif

