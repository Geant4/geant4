#ifndef _G4MAKEVOL_
#define _G4MAKEVOL_ 1

G4LogicalVolume* 
G4makevol(G4String& vname, G4String& shape, G4int nmed, G4double* Rpar, 
	  G4int npar);
#endif

