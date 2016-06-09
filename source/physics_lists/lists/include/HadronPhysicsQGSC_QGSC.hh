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
// $Id: HadronPhysicsQGSC_QGSC.hh,v 1.4 2009/04/14 07:23:08 mkossov Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_QGSC
//
// Author: 2009  M. Kosov based on HadronPhysicsQGSC_BERT
//
// Modified:
//
//----------------------------------------------------------------------------

#ifndef HadronPhysicsQGSC_QGSC_h
#define HadronPhysicsQGSC_QGSC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
//#include "G4MiscLHEPBuilder.hh"
#include "G4MiscQGSCBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4QGSC_QGSCPiKBuilder.hh"
//#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4QGSC_QGSCProtonBuilder.hh"
//#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4QGSC_QGSCNeutronBuilder.hh"
//#include "G4BertiniNeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"

class HadronPhysicsQGSC_QGSC : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSC_QGSC(const G4String& name ="hadron",G4bool quasiElastic=true);
    virtual ~HadronPhysicsQGSC_QGSC();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4QGSC_QGSCNeutronBuilder * theQGSCNeutron;
    //G4BertiniNeutronBuilder * theBertiniNeutron;
    G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission
    
    G4PiKBuilder * thePiK;
    G4QGSC_QGSCPiKBuilder * theQGSCPiK;
    //G4BertiniPiKBuilder * theBertiniPiK;
    
    G4ProtonBuilder * thePro;
    G4QGSC_QGSCProtonBuilder * theQGSCPro;    
    //G4BertiniProtonBuilder * theBertiniPro;

    G4MiscQGSCBuilder * theMiscQGSC;    
    //G4MiscLHEPBuilder * theMiscLHEP;
    
    G4bool QuasiElastic;
};

#endif

