// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gspos.cc,v 1.3 1999-05-12 08:10:07 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3RotTable.hh"

G4bool G3IsMany(G4String);

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

void G4gspos(G4String vname, G4int num, G4String vmoth, G4double x,
             G4double y, G4double z, G4int irot, G4String vonly)
{
  // get the logical volume pointer from the name
  G4LogicalVolume *lvol = G3Vol.GetLV(vname);
  
        // get the rotation matrix pointer from the G3 IROT index
  G3toG4RotationMatrix *rotm;
  if (irot>0) {
    rotm = G3Rot.get(irot);
  } else {
    rotm = NULL;
  }
  
  // translation offset
  G4ThreeVector *offset = new G4ThreeVector(x, y, z);
  
  // determine ONLY/MANY status
  G4bool isMany = G3IsMany(vonly);
  
  // check for negative parameters in volume definition
  G4double *pars = NULL;
  G4int npar;
  G4bool negpars = G3NegVolPars(pars,&npar,vname,vmoth,"GSPOS");
  
  // get the logical volume pointer of the mother from the name
  G4LogicalVolume *mothLV = G3Vol.GetLV(vmoth);
  G4PVPlacement* pvol = new G4PVPlacement((G4RotationMatrix*)rotm, 
					  *offset, lvol, vname,
					  mothLV, isMany, num);
  
  // add it to the List
  //    G3Vol.PutPV(&vname, pvol);
  delete offset;
}

