// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdvt.cc,v 1.3 1999-05-22 06:31:46 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3CalcParams.hh"

void PG4gsdvt(RWCString tokens[])
{
  // fill the parameter containers
  G3fillParams(tokens,PTgsdvt);
  
  // interpret the parameters
  G4String vname = Spar[0];
  G4String vmoth = Spar[1];
  G4int iaxis = Ipar[0];
  G4int numed = Ipar[1];
  G4int ndvmx = Ipar[2];
  G4double Step = Rpar[0];
  
  G4gsdvt(vname,vmoth,Step,iaxis,numed,ndvmx);
}

void G4gsdvt(G4String vname, G4String vmoth, G4double Step, G4int iaxis,
             G4int numed, G4int ndvmx){
  /*
  // get the physical volume pointer of the mother from the name
  G4int npv=0;
  G4VPhysicalVolume* mothPV;
  G3Vol.FindPV(&npv, &vmoth);
  if (npv == 0) {
    G4cout << " G4gsdvt error: No physical volume for " << vmoth << endl;
    return;
  }
  for (G4int ipv=0; ipv<npv; ipv++) {
    mothPV = G3Vol.GetPV(ipv);
    
    // extract the needed parameters from the mother's logical volume
    G4double rangehi;
    G4double rangelo;
    EAxis axiscode;
    G4String shape;
    G4int nmed;
    G4double *Rpar = NULL;
    G4int npar;
    G4int zeronpar=0;
    G4VSolid *solid = NULL;
    G3Vol.GetLVPars(&vmoth, iaxis, &rangehi, &rangelo, &axiscode, &shape,
                    &nmed, &Rpar, &npar, solid);
    // Calculate the number of divisions
    G4double c0 = rangelo;
    G4int ndiv = int((rangehi - rangelo) / Step);
    if (ndvmx > 255) ndvmx = 255;
    if (ndiv > ndvmx && ndvmx > 0 ) {
      ndiv = ndvmx;
      c0 = 0.5*((rangehi - rangelo) - ndiv*Step) + rangelo;
    }
    // nullify parameter values
    for (G4int i=0; i<npar; i++) Rpar[i] = 0.;
    // Generate the logical volume for the subvolumes
    G4gsvolu(vname, shape, numed, Rpar, zeronpar);
    // and obtain the pointer
    G4LogicalVolume *lvol = G3Vol.GetLVx(vname);
    
    // check for negative parameters in volume definition
    G4double *pars = NULL;
    //   G4bool negpars = G3NegVolPars(pars,&npar,vname,vmoth,"GSDVN");
    G4bool negpars = 0;
    G4double width = rangehi - c0;
    G4double offset = (rangehi + c0)/2.;
    
    if ( ! negpars ) {
      // Generate replicas
      G4PVReplica *pvol = new G4PVReplica(
         vname, lvol, mothPV, axiscode, ndiv, width, offset);
      G3Vol.PutPV1(&vname, pvol);
    } else {
      G4VPVParameterisation *dvnParam = new G3CalcParams(ndiv, width, offset);
      G4PVParameterised *pvol = new G4PVParameterised(
         vname, lvol, mothPV, axiscode, ndiv, dvnParam);
      G3Vol.PutPV1(&vname, pvol);
    }
  }
  */
  G4cerr << "G4gsdvt currently not supported" << endl;
}
