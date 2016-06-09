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
// $Id: HadronPhysicsLHEP_GN.hh,v 1.4 2006/06/29 17:56:47 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsLHEP_GN
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: don't  keep processes as data members, but new these
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsLHEP_GN_h
#define HadronPhysicsLHEP_GN_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4LHEPPiKBuilder.hh"
#include "G4PiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LHEPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LHEPNeutronBuilder.hh"

class HadronPhysicsLHEP_GN : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsLHEP_GN(const G4String& name ="hadron");
    virtual ~HadronPhysicsLHEP_GN();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4LHEPNeutronBuilder * theLHEPNeutron;
    
    G4PiKBuilder * thePiK;
    G4LHEPPiKBuilder * theLHEPPiK;
    
    G4ProtonBuilder * thePro;
    G4LHEPProtonBuilder * theLHEPPro;
    
    G4MiscLHEPBuilder * theMiscLHEP;
};

#endif

