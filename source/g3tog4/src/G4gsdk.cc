//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4gsdk.cc,v 1.5 2001-07-11 09:59:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4Decay.hh"
#include "G3toG4.hh"
#include "G3PartTable.hh"

void PG4gsdk(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdk);

    // interpret the parameters
    G4int ipart = Ipar[0];
    G4int *mode = &Ipar[3];
    G4double *bratio = Rpar;

    G4gsdk(ipart,bratio,mode);
}

void G4gsdk(G4int ipart, G4double bratio[], G4int mode[])
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
