// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsposp.cc,v 1.1 1999-01-07 16:06:51 gunter Exp $
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
G4LogicalVolume* G4makevol(G4String vname, G4String shape, G4int nmed,
                           G4double Rpar[], G4int npar);

void PG4gsposp(RWCString tokens[])
{
        // fill the parameter containers
    G3fillParams(tokens,PTgsposp);
  
  // interpret the parameters
    G4String name = Spar[0];
    G4String moth = Spar[1];
    G4String only = Spar[2];
    G4int num = Ipar[0];
    G4int irot = Ipar[1];
    G4int npar = Ipar[2];
    G4double x = Rpar[0];
    G4double y = Rpar[1];
    G4double z = Rpar[2];
    G4double *pars = &Rpar[3];
  
    G4gsposp(name, num, moth, x, y, z, irot, only, pars, npar);
}

void G4gsposp(G4String vname, G4int num, G4String vmoth, G4double x,
              G4double y, G4double z, G4int irot, G4String vonly,
              G4double pars[], G4int npar)
{
    G4bool _debug=false;
    
        // get the rotation matrix pointer from the G3 IROT index
    G4RotationMatrix* rotm = G3Rot.get(irot);
  
        // translation offset
    G4ThreeVector* offset = new G4ThreeVector(x, y, z);
  
        // determine ONLY/MANY status
    G4bool isMany = G3IsMany(vonly);
  
        // check for negative parameters in volume definition.
        // pars is updated if negative params are found.
    G4bool negpars = G3NegVolPars(pars,&npar,vname,vmoth,"GSPOS");

        // get the logical volume pointer of the mother from the name
    G4LogicalVolume *mothLV = G3Vol.GetLVx(vmoth);
    
        // create a logical volume and add it to the constituent List of the
        // G3 logical volume
    G4String shape;
    G4int nmed;
    G3Vol.GetLVInfo(&vname, &shape, &nmed);
    if (_debug) {
        G4cout << "gsposp: creating lvol " << vname << " shape " << shape << " nmed " <<
            nmed << " moth " << vmoth << " mothLV " << mothLV << " pars ";
        for (int i=0; i<npar; i++) G4cout << pars[i] << " ";
        G4cout << endl;
    }      
    G4LogicalVolume* lvol = G4makevol(vname, shape, nmed, pars, npar);
        // add the logical volume to the constituent List
    G3Vol.AddConstituentLVol(&vname, lvol);
    
    G4VPhysicalVolume* pvol = new G4PVPlacement(rotm, *offset, lvol, vname,
                                                mothLV, isMany, num);

        // add it to the List
    G3Vol.PutPV(&vname, pvol);
    delete offset;
}
