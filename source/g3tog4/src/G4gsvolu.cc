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
// $Id: G4gsvolu.cc,v 1.8 2001-07-11 09:59:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 13.10.99

#include "g4std/iomanip"
#include "G3VolTable.hh"
#include "globals.hh"
#include "G3toG4.hh"
#include "G3toG4MakeSolid.hh"

void PG4gsvolu(G4String tokens[]) {
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
                             G4double Rpar[], G4int npar)
{    
  // create the solid
  G4bool hasNegPars;
  G4bool deferred;   
  G4bool okAxis[3];
  G4VSolid* solid
    = G3toG4MakeSolid(vname, shape, Rpar, npar, hasNegPars, deferred, okAxis);  

  // if solid has been deferred 
  // VTE is created with hasNegPars = true  
  if (deferred) hasNegPars = true;   

  // create VTE
  G3VolTableEntry* vte 
     = new G3VolTableEntry(vname, shape, Rpar, npar, nmed, solid, hasNegPars);
  G3Vol.PutVTE(vte);
  
  return vte;
}

void G4gsvolu(G4String vname, G4String shape, G4int nmed, G4double* Rpar,
              G4int npar)
{
  /*
  G4cout << "Creating logical volume " << vname << " shape " << shape
  	 << " nmed " << nmed << " #pars "<< npar << " parameters (cm): ";
  for (int ipar=0; ipar< npar; ipar++) G4cout << G4std::setw(8) << Rpar[ipar];
  G4cout << G4endl;
  */
  if (G3Vol.GetVTE(vname)) {
    // abort if VTE with given name exists
    G4Exception("G4gsvolu: Attempt to create volume " + vname + " twice.");
  }
  else {  
    G4CreateVTE(vname, shape, nmed, Rpar, npar);
  }  
}
