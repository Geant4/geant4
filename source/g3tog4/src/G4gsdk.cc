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
// $Id: G4gsdk.cc 67982 2013-03-13 10:36:03Z gcosmo $

#include "G4Decay.hh"
#include "G3toG4.hh"
#include "G3PartTable.hh"

void PG4gsdk(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdk);

    // interpret the parameters
    G4int ipart = Ipar[0];
    G4int *mode = &Ipar[3];
    G4double *bratio = Rpar;

    G4gsdk(ipart,bratio,mode);
}

void G4gsdk(G4int, G4double*, G4int*)
{
/*
    // create decay object for the particle
    G4Decay *decay = new G4Decay();
    // add decay modes
    for (G4int i=0; i<6; i++) {
        if (mode[i] != 0) {
// $$$            decay->AddMode(mode[i],bratio[i]);
        }
    }
    // associate decay object with particle ipart
    G4ParticleDefinition *part = G3Part.Get(ipart);
// $$$    part->SetDecay(decay);
*/
}
