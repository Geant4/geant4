// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdvn2.cc,v 1.5 1999-11-15 10:39:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3CalcParams.hh"

void PG4gsdvn2(G4String tokens[]) {
  // fill the parameter containers
  G3fillParams(tokens, PTgsdvn2);
  
  // interpret the parameters
  G4String vname = Spar[0];
  G4String vmoth = Spar[1];
  G4int ndiv = Ipar[0];
  G4int iaxis = Ipar[1];
  G4int numed = Ipar[2];
  G4double c0 = Rpar[0];
  
  G4gsdvn2(vname, vmoth, ndiv, iaxis, c0, numed);
}

void G4gsdvn2(G4String vname, G4String vmoth, G4int ndiv, G4int iaxis,
              G4double c0, G4int numed){
  /*
  // get the physical volume pointer of the mother from the name
  G4int npv=0;
  G4VPhysicalVolume* mothPV;
  G3Vol.FindPV(&npv, &vmoth);
  if (npv == 0) {
    G4cout << " G4gsdvn2 error: No physical volume for " << vmoth << endl;
    return;
  }
  for (G4int ipv=0; ipv<npv; ipv++) {
    mothPV = G3Vol.GetPV(ipv);
  
    // extract the needed parameters from the mother's logical volume
    G4double rangehi=0;
    G4double rangelo=0;
    G4String shape;
    G4int nmed;
    G4double *Rpar = NULL;
    G4int npar;
    G4int zeronpar=0;
    G4VSolid *solid = NULL;
    G3Vol.GetLVPars(&vmoth, &shape, &nmed, &Rpar, &npar, solid);
    // nullify parameter values
    for (G4int i=0; i<npar; i++) Rpar[i] = 0.;
    // Generate the logical volume for the subvolumes
    G4gsvolu(vname, shape, numed, Rpar, zeronpar);
    // and obtain the pointer
    G4LogicalVolume *lvol = G3Vol.GetLVx(vname);
    // check for negative parameters in volume definition
    G4double *pars = 0;
    //    G4bool negpars = G3NegVolPars(pars, &npar, vname, vmoth, "GSDVN");
    G4bool negpars = false;
    G4double width = rangehi - c0;
    G4double offset = (rangehi + c0)/2.;
  
    if (! negpars ) {
      Generate replicas
        G4PVReplica *pvol = new 
	G4PVReplica(vname, lvol, mothPV, axiscode, ndiv, width, offset);
      G3Vol.PutPV1(&vname, pvol);
    } else {
      G4VPVParameterisation *dvnParam = new G3CalcParams(ndiv, width, offset);
      G4PVParameterised *pvol = new 
	G4PVParameterised(vname, lvol, mothPV, axiscode, ndiv, dvnParam);
      G3Vol.PutPV1(&vname, pvol);
    }
  }
  */
  G4cerr << "G4gsdvn2 currently not supported." << endl;
}

  
