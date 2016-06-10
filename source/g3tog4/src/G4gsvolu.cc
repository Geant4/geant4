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
// $Id: G4gsvolu.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 13.10.99

#include <iomanip>
#include "G3VolTable.hh"
#include "globals.hh"
#include "G3toG4.hh"
#include "G3toG4MakeSolid.hh"

void PG4gsvolu(G4String *tokens) {
    // fill the parameter containers
    G3fillParams(tokens,PTgsvolu);

    // interpret the parameters
    G4String vname = Spar[0];
    G4String shape = Spar[1];
    G4int nmed = Ipar[0];
    G4int npar = Ipar[1];
    G4double *pars = Rpar;

    G4gsvolu(vname, shape, nmed, pars, npar);
}

G3VolTableEntry* G4CreateVTE(G4String vname, G4String shape, G4int nmed,
                             G4double rpar[], G4int npar)
{    
  // create the solid
  G4bool hasNegPars;
  G4bool deferred;   
  G4bool okAxis[3];
  G4VSolid* solid
    = G3toG4MakeSolid(vname, shape, rpar, npar, hasNegPars, deferred, okAxis);  

  // if solid has been deferred 
  // VTE is created with hasNegPars = true  
  if (deferred) hasNegPars = true;   

  // create VTE
  G3VolTableEntry* vte 
     = new G3VolTableEntry(vname, shape, rpar, npar, nmed, solid, hasNegPars);
  G3Vol.PutVTE(vte);
  
  return vte;
}

void G4gsvolu(G4String vname, G4String shape, G4int nmed, G4double* rpar,
              G4int npar)
{
  /*
  G4cout << "Creating logical volume " << vname << " shape " << shape
  	 << " nmed " << nmed << " #pars "<< npar << " parameters (cm): ";
  for (int ipar=0; ipar< npar; ipar++) G4cout << std::setw(8) << rpar[ipar];
  G4cout << G4endl;
  */
  if (G3Vol.GetVTE(vname)) {
    // abort if VTE with given name exists
    G4String text = "G4gsvolu: Attempt to create volume " + vname + " twice.";
    G4Exception("G4gsvolu()", "G3toG40024", FatalException, text);
    return;
  }
  else {  
    G4CreateVTE(vname, shape, nmed, rpar, npar);
  }  
}
