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
// $Id: G4gsdeth.cc 67982 2013-03-13 10:36:03Z gcosmo $

#include "G3toG4.hh"
#include "G3DetTable.hh"

class G4VSensitiveDetector;

void PG4gsdeth(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdeth);

    // interpret the parameters
    G4String chset = Spar[0];
    G4String chdet = Spar[1];
    G4int nh = Ipar[0];
    G4String chnamh[100];
    for (G4int i=0; i<=nh; i++ ) chnamh[i] = Spar[2+i].data();
    G4int *nbitsh = &Ipar[1];
    G4double *orig = Rpar;
    G4double *fact = &Rpar[nh];

    G4gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact);
}

void G4gsdeth(G4String, G4String, G4int, G4String*,
              G4int*, G4double*, G4double*)
{
    // Get pointer to sensitive detector chset
    // G4VSensitiveDetector* sdet = G3Det.GetSD(chset, chdet);
    // Add hits to sensitive detector
    // for (G4int i=0; i<nh; i++) {
      // $$$        sdet->AddHit(chnamh[i],nbitsh[i],orig[i],fact[i]);
    // }
}
