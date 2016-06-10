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
// $Id: G4gsdvx.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, V.Berejnoi, 27 Sep 99

#include "globals.hh"
#include "G3toG4.hh"

void G4gsdvn2(G4String name, G4String moth, G4int ndiv, G4int iaxis,
              G4double c0, G4int numed);

void G4gsdvt2(G4String name, G4String moth, G4double Step, G4int iaxis,
              G4double c0, G4int numed, G4int ndvmx);

void PG4gsdvx(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdvx);

    // interpret the parameters
    G4String name = Spar[0];
    G4String moth = Spar[1];
    G4int ndiv = Ipar[0];
    G4int iaxis = Ipar[1];
    G4int numed = Ipar[2];
    G4int ndvmx = Ipar[3];
    G4double Step = Rpar[0];
    G4double c0 = Rpar[1];

    G4gsdvx(name,moth,ndiv,iaxis,Step,c0,numed,ndvmx);
}

void G4gsdvx(G4String name, G4String moth, G4int ndiv, G4int iaxis,
             G4double Step, G4double c0, G4int numed, G4int ndvmx)
{
    // pass to gsdvn2 or gsdvt2
    if (Step > 0.) {
        G4gsdvt2(name,moth,Step,iaxis,c0,numed,ndvmx);
    } else if (ndiv > 0) {
        G4gsdvn2(name,moth,ndiv,iaxis,c0,numed);
    }
}
