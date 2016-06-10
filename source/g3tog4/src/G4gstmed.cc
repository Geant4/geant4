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
// $Id: G4gstmed.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// The last G4int argument of G4gstmed(..) is used for sending
// info whether the Geant3 tracking medium parameters should
// be set (the max step, later: G3 default cut values).
//
// by I.Hrivnacova, 27 Sep 99

#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"
#include "G4UserLimits.hh"
#include "G4MagneticField.hh"
#include "G4Material.hh"

void PG4gstmed(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgstmed);

    // interpret the parameters
    G4String name = Spar[0];
    G4int itmed = Ipar[0];
    G4int nmat = Ipar[1];
    G4int isvol = Ipar[2];
    G4int ifield = Ipar[3];
    G4int nwbuf = Ipar[4];
    G4double fieldm = Rpar[0];
    G4double tmaxfd = Rpar[1];
    G4double stemax = Rpar[2];
    G4double deemax = Rpar[3];
    G4double epsil = Rpar[4];
    G4double stmin = Rpar[5];
    G4double *ubuf = &Rpar[6];

    G4gstmed(itmed,name,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
             deemax,epsil,stmin,ubuf,nwbuf);
}

void G4gstmed(G4int itmed, G4String, G4int nmat, G4int isvol,
              G4int, G4double, G4double,
              G4double stemax, G4double, G4double,
              G4double, G4double*, G4int useG3TMLimits)
{
    // get the pointer to material nmat
    G4Material* material = G3Mat.get(nmat);

    // NB. there is the possibility for redundancy in the mag field
    //     and user limits objects. Who cares.
    // Generate the mag field object
    // $$$ G4MagneticField* field = new G4MagneticField(ifield, fieldm, tmaxfd);
    G4MagneticField* field = 0;

    // Generate the user limits object
    // !!! only "stemax" has its equivalent in G4

    G4UserLimits* limits = 0;
    if (useG3TMLimits) {
      limits = new G4UserLimits();
      limits->SetMaxAllowedStep(stemax*cm);
      // limits->SetG3DefaultCuts();
         // this is arranged globally by physics manager
    }

    // Store this medium in the G3Med structure
    G3Med.put(itmed, material, field, limits, isvol);
}
