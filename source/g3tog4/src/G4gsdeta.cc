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
// $Id: G4gsdeta.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
//
// modified by I.Hrivnacova, 6 Oct 99

#include "G3toG4.hh"
#include "G3DetTable.hh"

void G4gsdeta(G4String chset, G4String chdet, G4String,
              G4int nwhi, G4int nwdi);

void PG4gsdeta(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdeta);

    // interpret the parameters
    G4String chset = Spar[0];
    G4String chdet = Spar[1];
    G4String chali = Spar[2];
    G4int nwhi = Ipar[0];
    G4int nwdi = Ipar[1];

    G4gsdeta(chset,chdet,chali,nwhi,nwdi);
}

void G4gsdeta(G4String chset, G4String chdet, G4String,
              G4int nwhi, G4int nwdi)
{
    G4int idtyp = G3Det.GetID(chset, chdet);
    // just associate another sensitive detector structure with
    // the volume chdet
    G4gsdetv(chset, chdet, idtyp, nwhi, nwdi);
}
