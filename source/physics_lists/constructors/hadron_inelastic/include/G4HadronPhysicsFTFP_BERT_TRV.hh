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
// $Id: G4HadronPhysicsFTFP_BERT_TRV.hh 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007  Gunter Folger
//   created from G4HadronPhysicsFTFP
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
// 01.11.2012 W.Pokorski & A.Ribon: use new cross sections
// 18.07.2016 A.Dotti: refactor following new code
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTFP_BERT_TRV_h
#define G4HadronPhysicsFTFP_BERT_TRV_h 1

#include "G4HadronPhysicsFTFP_BERT.hh"

class G4HadronPhysicsFTFP_BERT_TRV : public G4HadronPhysicsFTFP_BERT
{
  public: 
    G4HadronPhysicsFTFP_BERT_TRV(G4int verbose =1);
    G4HadronPhysicsFTFP_BERT_TRV(const G4String& name, G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsFTFP_BERT_TRV() {}

  private:
    //Modify the minimum needed
      virtual void Pion() override;
      virtual void Kaon() override;
      virtual void DumpBanner() override;
};

#endif
