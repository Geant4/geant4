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
// $Id: G4HadronPhysicsQGSP_FTFP_BERT.hh 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_FTFP_BERT
//
// Author: 2 Apr 2009  J.Apostolakis & V.Ivantchenko
//
// Replace the use of LEP for proton, neutron, pions
//
// Modified:
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGSP_FTFP_BERT_h
#define G4HadronPhysicsQGSP_FTFP_BERT_h 1


#include "globals.hh"
#include "G4ios.hh"

#include "G4HadronPhysicsQGSP_BERT.hh"


class G4HadronPhysicsQGSP_FTFP_BERT : public G4HadronPhysicsQGSP_BERT 
{
  public: 
    G4HadronPhysicsQGSP_FTFP_BERT(G4int verbose =1);
    G4HadronPhysicsQGSP_FTFP_BERT(const G4String& name, G4bool quasiElastic=true);
  protected:
    virtual void DumpBanner() override;    
    G4bool QuasiElastic;
};


#endif


