// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdetv.cc,v 1.1 1999-01-07 16:06:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3toG4.hh"
#include "G3DetTable.hh"
#include "G3VolTable.hh"

class G4VSensitiveDetector;

void PG4gsdetv(RWCString tokens[])
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

void G4gsdetv(G4String chset, G4String chdet, G4int idtyp, G4int,
              G4int)
{
    // get lvol for detector chdet
    G4LogicalVolume *lvol = G3Vol.GetLVx(chdet);
    if (lvol == NULL) {
      G4cout << "G4gsdetv: Logical volume " << chdet << " not available. Skip." << endl;
      return;
    }
    G4cout << "G4gsdetv not currently implemented." << endl;
    
    // Generate a sensitive detector structure
        // G4VSensitiveDetector *sdet;
    // $$$    G4VSensitiveDetector *sdet = new G4VSensitiveDetector(chset);
    // inform the logical volume of its sensitive detector
        // lvol->SetSensitiveDetector(sdet);
    // $$$ sdet->SetID(idtyp);
    // Add the sensitive detector to the table
        // G3Det.put(chset,idtyp,sdet);
}

