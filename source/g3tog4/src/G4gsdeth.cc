// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdeth.cc,v 1.1 1999-01-07 16:06:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3toG4.hh"
#include "G3DetTable.hh"

class G4VSensitiveDetector;

void PG4gsdeth(RWCString tokens[])
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

void G4gsdeth(G4String chset, G4String, G4int nh, G4String chnamh[],
              G4int nbitsh[], G4double orig[], G4double fact[])
{
    // Get pointer to sensitive detector chset
    G4VSensitiveDetector* sdet = G3Det.get(chset);
    // Add hits to sensitive detector
    for (G4int i=0; i<nh; i++) {
      // $$$        sdet->AddHit(chnamh[i],nbitsh[i],orig[i],fact[i]);
    }
}
