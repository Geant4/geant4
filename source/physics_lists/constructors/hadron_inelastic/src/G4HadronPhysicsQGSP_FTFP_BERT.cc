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
// ClassName:   G4HadronPhysicsQGSP_FTFP_BERT
//
// Authors: 2 Apr 2009 J.Apostolakis/V.Ivantchenko: created starting from QGSP_BERT
//
// Modified:
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_FTFP_BERT);


G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(G4int verb)
    : G4HadronPhysicsQGSP_FTFP_BERT("hInelastic QGSP_FTFP_BERT",true) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(const G4String& name, 
							 G4bool quasiElastic)
    : G4HadronPhysicsQGSP_BERT(name,quasiElastic)
{
  G4HadronicParameters::Instance()->SetEnableBCParticles(true);
}

G4HadronPhysicsQGSP_FTFP_BERT::~G4HadronPhysicsQGSP_FTFP_BERT()
{}

