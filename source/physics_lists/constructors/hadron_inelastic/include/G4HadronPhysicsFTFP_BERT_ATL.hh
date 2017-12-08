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
// $Id:$
//
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   April 2016
//
// Hadron physics for the new physics list FTFP_BERT_ATL.
// This is a modified version of the FTFP_BERT hadron physics for ATLAS.
// The hadron physics of FTFP_BERT_ATL has the transition between Bertini
// (BERT) intra-nuclear cascade model and Fritiof (FTF) string model in the
// energy region [9, 12] GeV (instead of [4, 5] GeV as in FTFP_BERT).
//
// Modified:
// 18.07.2017 A.Dotti: refactor forllowin new code
//---------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTFP_BERT_ATL_h
#define G4HadronPhysicsFTFP_BERT_ATL_h 1

#include "G4HadronPhysicsFTFP_BERT.hh"

class G4HadronPhysicsFTFP_BERT_ATL : public G4HadronPhysicsFTFP_BERT
{
  public: 
    G4HadronPhysicsFTFP_BERT_ATL(G4int verbose =1);
    G4HadronPhysicsFTFP_BERT_ATL(const G4String& name, G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsFTFP_BERT_ATL() {}

  private:
    //Modify the minimum needed
    virtual void Pion() override;
    virtual void Kaon() override;
    virtual void DumpBanner() override;
};

#endif

