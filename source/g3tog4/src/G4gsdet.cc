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
// $Id: G4gsdet.cc 67982 2013-03-13 10:36:03Z gcosmo $

#include "globals.hh"
#include "G3toG4.hh"

void PG4gsdet(G4String* tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdet);

    // interpret the parameters
    G4String chset = Spar[0];
    G4String chdet = Spar[1];
    G4int nv = Ipar[0];
    G4String chnmsv[100];
    for (G4int i=0; i<=nv; i++ ) chnmsv[i] = Spar[2+i].data();
    G4int *nbits = &Ipar[1];
    G4int idtyp = Ipar[1+nv];
    G4int nwhi = Ipar[2+nv];
    G4int nwdi = Ipar[3+nv];

    G4gsdet(chset,chdet,nv,chnmsv,nbits,idtyp,nwhi,nwdi);
}

void G4gsdet(G4String chset, G4String chdet, G4int, G4String*,
             G4int*, G4int idtyp, G4int nwhi, G4int nwdi)
{ 
    G4gsdetv(chset, chdet, idtyp, nwhi, nwdi);
}
