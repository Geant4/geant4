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
// $Id: G4gspart.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
#include "G4ProcessManager.hh"
#include "G3toG4.hh"
#include "G3PartTable.hh"

void PG4gspart(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgspart);

    // interpret the parameters
    G4String chnpar = Spar[0];
    G4int ipart = Ipar[0];
    G4int itrtyp = Ipar[1];
    G4int nwb = Ipar[2];
    G4double amass = Rpar[0];
    G4double charge = Rpar[1];
    G4double tlife = Rpar[2];
    G4double *ubuf = &Rpar[3];

    G4gspart(ipart,chnpar,itrtyp,amass,charge,tlife,ubuf,nwb);
}

void G4gspart(G4int, G4String, G4int, G4double,
              G4double, G4double, G4double*, G4int)
{
}

#ifdef OMIT_CODE
void G4gspart(G4int ipart, G4String chnpar, G4int itrtyp, G4double amass,
              G4double charge, G4double tlife, G4double*, G4int)
{
    // Handle conversion of itrtyp into an appropriate ProcessManager
    G4ProcessManager *mgr = 0;
    switch (itrtyp) {
    case 1:
        // gamma
// $$$        mgr = gammaProcessManager;
        break;
    case 2:
        // electron
// $$$        mgr = electronProcessManager;
        break;
    case 3:
        // neutron
// $$$        mgr = neutronProcessManager;
        break;
    case 4:
        // hadron
// $$$        mgr = hadronProcessManager;
        break;
    case 5:
        // muon
// $$$        mgr = muonProcessManager;
        break;
    case 6:
        // geantino
// $$$        mgr = geantinoProcessManager;
        break;
    case 7:
        // heavy ion
// $$$        mgr = heavyIonProcessManager;
        break;
    case 8:
        // light ion
// $$$        mgr = lightIonProcessManager;
        break;
    default:
// $$$        mgr = theProcessManager;
        break;
    }

    // Create the particle; provide parameters and a process mgr.
    G4ParticleDefinition *part = new G4ParticleDefinition(chnpar);
    part->SetProcessManager(mgr);
    part->SetPDGMass(amass);
    part->SetPDGCharge(charge);
    part->SetPDGLifeTime(tlife);

    // add to particle table
    G3Part.put(&ipart, part);
}
#endif
