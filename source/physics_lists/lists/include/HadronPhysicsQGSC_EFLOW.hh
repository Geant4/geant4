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
// $Id: HadronPhysicsQGSC_EFLOW.hh,v 1.2 2007/04/26 14:47:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_EFLOW
//
// Author: 2006 Gunter Folger
//
// Created from HadronPhysicsQGSC
// Modified:
// 25.04.2007 G.Folger: Add quasielastic as option, use quasielastic by default
//
//----------------------------------------------------------------------------

#ifndef HadronPhysicsQGSC_EFLOW_h
#define HadronPhysicsQGSC_EFLOW_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4LEPPiKBuilder.hh"
#include "G4QGSCEflowPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LEPProtonBuilder.hh"
#include "G4QGSCEflowProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"
#include "G4QGSCEflowNeutronBuilder.hh"

class HadronPhysicsQGSC_EFLOW : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSC_EFLOW(const G4String& name ="hadron",G4bool quasiElastic=true);
    virtual ~HadronPhysicsQGSC_EFLOW();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4LEPNeutronBuilder * theLEPNeutron;
    G4QGSCEflowNeutronBuilder * theQGSCEflowNeutron;
    
    G4PiKBuilder * thePiK;
    G4LEPPiKBuilder * theLEPPiK;
    G4QGSCEflowPiKBuilder * theQGSCEflowPiK;
    
    G4ProtonBuilder * thePro;
    G4LEPProtonBuilder * theLEPPro;
    G4QGSCEflowProtonBuilder * theQGSCEflowPro;    
    
    G4MiscLHEPBuilder * theMiscLHEP;
    
    G4bool QuasiElastic;
};

#endif

