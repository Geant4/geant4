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
// $Id: G4gsdetu.cc,v 1.5 2003/04/03 11:30:55 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $

#include "G3toG4.hh"

void PG4gsdetu(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdetu);

    // interpret the parameters
    G4String chset = Spar[0].data();
    G4String chdet = Spar[1].data();
    G4int nupar = Ipar[0];
    G4double *upar = Rpar;

    G4gsdetu(chset,chdet,nupar,upar);
}

void G4gsdetu(G4String, G4String, G4int, G4double*)
{
    // $$$ nothing right now
}
