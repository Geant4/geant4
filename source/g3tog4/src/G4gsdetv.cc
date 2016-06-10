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
// $Id: G4gsdetv.cc 67982 2013-03-13 10:36:03Z gcosmo $

#include "G4ios.hh"
#include "G3toG4.hh"
#include "G3DetTable.hh"
#include "G3VolTable.hh"

class G4VSensitiveDetector;

void PG4gsdetv(G4String *tokens)
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdetv);

    // interpret the parameters
    G4String chset = Spar[0];
    G4String chdet = Spar[1];
    G4int idtyp = Ipar[0];
    G4int nwhi = Ipar[1];
    G4int nwdi = Ipar[2];

    G4gsdetv(chset,chdet,idtyp,nwhi,nwdi);
}

void G4gsdetv(G4String, G4String, G4int, G4int,G4int)
{  
  G4cout << "G4gsdetv not currently implemented." << G4endl;
  /*
    // get lvol for detector chdet
    G4LogicalVolume *lvol = G3Vol.GetLV(chdet);
    if (lvol == 0) {
    G4cout << "G4gsdetv: Logical volume " << chdet << " not available. Skip." << G4endl;
    return;
    }
    // Generate a sensitive detector structure
    // G4VSensitiveDetector *sdet;
    // $$$    G4VSensitiveDetector *sdet = new G4VSensitiveDetector(chset);
    // inform the logical volume of its sensitive detector
    // lvol->SetSensitiveDetector(sdet);
    // $$$ sdet->SetID(idtyp);
    // Add the sensitive detector to the table
    // G3Det.put(chset,idtyp,sdet);
    */
}

