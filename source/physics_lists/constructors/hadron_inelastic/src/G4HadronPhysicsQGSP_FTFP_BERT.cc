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
// $Id: G4HadronPhysicsQGSP_FTFP_BERT.cc 105736 2017-08-16 13:01:11Z gcosmo $
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
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_FTFP_BERT);


G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(G4int)
    : G4HadronPhysicsQGSP_FTFP_BERT("hInelastic QGSP_FTFP_BERT",true) {}

G4HadronPhysicsQGSP_FTFP_BERT::G4HadronPhysicsQGSP_FTFP_BERT(const G4String& name, 
							 G4bool quasiElastic)
    :  G4HadronPhysicsQGSP_BERT(name,quasiElastic), QuasiElastic(quasiElastic)
{
 maxBERT_proton = maxBERT_neutron = maxBERT_pik = 8.*GeV;
 minFTFP_proton = minFTFP_neutron = minFTFP_pik = 6.*GeV;
}


void G4HadronPhysicsQGSP_FTFP_BERT::DumpBanner()
{
  G4cout << " New QGSP_FTFP_BERT physics list, replaces LEP with FTF/P for p/n/pi (/K?)";
  G4cout << "  Thresholds: " << G4endl;
  G4cout << "    1) between BERT  and FTF/P over the interval " 
	 << minFTFP_proton/GeV << " to " << maxBERT_proton/GeV << " GeV. " << G4endl;
  G4cout << "    2) between FTF/P and QGS/P over the interval " 
	 << minQGSP_proton/GeV << " to " << maxFTFP_proton/GeV << " GeV. " << G4endl;
  G4cout << "  -- quasiElastic was asked to be " << QuasiElastic << G4endl
	 << "     Changed to " << QuasiElasticQGS << " for QGS "
	 << " and to " << QuasiElasticFTF << " (must be false) for FTF" << G4endl;
}
