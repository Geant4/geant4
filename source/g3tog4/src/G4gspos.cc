// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gspos.cc,v 1.7 1999-05-22 06:31:52 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "VolTableEntry.hh"
#include "G3Pos.hh"

void PG4gspos(RWCString tokens[])
{
        // fill the parameter containers
    G3fillParams(tokens,PTgspos);
  
        // interpret the parameters
    G4String name = Spar[0];
    G4String moth = Spar[1];
    G4String only = Spar[2];
    G4int num = Ipar[0];
    G4int irot = Ipar[1];
    G4double x = Rpar[0]*cm;
    G4double y = Rpar[1]*cm;
    G4double z = Rpar[2]*cm;
  
    G4gspos(name, num, moth, x, y, z, irot, only);
}

void G4gspos(const G4String& vname, const G4int num, const G4String& vmoth, 
	     G4double x, G4double y, G4double z, const G4int irot, 
	     const G4String& vonly){

  VolTableEntry* _VTE = G3Vol.GetVTE(vname);
  VolTableEntry* MVTE = G3Vol.GetVTE(vmoth);

  if (_VTE == 0) {
    G4cerr << "G4gspos: '" << vname << "' has no VolTableEntry" << endl;
  } else if (MVTE == 0) {
    G4cerr << "G4gspos: '" << vname << "' mother volume '" << vmoth 
	   << "' has no VolTableEntry" << endl;
  } else {
    
    // translation offset
    G4ThreeVector* offset = new G4ThreeVector(x, y, z);
    
    // create a G3Pos object and add it to the G3VolTable

    G3Pos* aG3Pos = new G3Pos(vname, num, vmoth, offset, irot, vonly);
    _VTE->AddG3Pos(aG3Pos);

    // add the G3Pos object as a daughter to its mother
    MVTE->AddDaughter(aG3Pos);
  }
}













