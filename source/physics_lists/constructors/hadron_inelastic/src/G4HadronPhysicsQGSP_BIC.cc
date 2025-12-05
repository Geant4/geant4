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
// ClassName:   G4HadronPhysicsQGSP_BIC
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
// 05.05.2020 A.Ribon: Use QGSP for antibaryons at high energies
// 07.05.2020 A.Ribon: Use QGSP for hyperons (and anti-hyperons) at high energies
//
//----------------------------------------------------------------------------
//

#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4HadronicParameters.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC);

G4HadronPhysicsQGSP_BIC::G4HadronPhysicsQGSP_BIC(G4int verb)
  : G4HadronPhysicsQGSP_BIC("hInelastic QGSP_BIC", true) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BIC::G4HadronPhysicsQGSP_BIC(const G4String& name, G4bool b)
  : G4HadronPhysicsQGSP_BERT(name, b)
{
  maxBIC_proton = maxBIC_neutron = 1.5*CLHEP::GeV;
  minBERT_proton = minBERT_neutron = 1.0*CLHEP::GeV;
  G4HadronicParameters::Instance()->SetEnableBCParticles(false);
  G4HadronicParameters::Instance()->SetUseRFilesForXS(true);
}

