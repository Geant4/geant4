// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3G4Interface.hh,v 1.6 1999-12-09 00:04:58 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//   Interfaces for G3 equivalent routines
//

#include "globals.hh"

class G4LogicalVolume;

void G4gsvolu(G4String name, G4String shape, G4int nmed, G4double* par,
              G4int npar);

void G4gspos(G4String name, G4int num, G4String moth, 
	     G4double x, G4double y, G4double z, G4int irot, 
	     G4String only);

void G4gsposp(G4String name, G4int num, G4String moth, 
	      G4double x, G4double y, G4double z, G4int irot, 
	      G4String only, G4double Rpar[], G4int npar);

void G4gsrotm(G4int irot, G4double theta1, G4double phi1,
              G4double theta2, G4double phi2, G4double theta3, G4double phi3);

void G4gsatt(G4String name, G4String attr, G4int ival);

void G4gsdvn(G4String vname, G4String vmoth, G4int ndiv, G4int iaxis);

void G4gsdvt(G4String name, G4String moth, G4double Step, G4int iaxis,
             G4int numed, G4int ndvmx);

void G4gsdvx(G4String name, G4String moth, G4int ndiv, G4int iaxis,
             G4double Step, G4double c0, G4int numed, G4int ndvmx);

void G4gsdvn2(G4String name, G4String moth, G4int ndiv, G4int iaxis,
              G4double c0, G4int numed);

void G4gsdvt2(G4String name, G4String moth, G4double Step, G4int iaxis,
              G4double c0, G4int numed, G4int ndvmx);

void G4gsmate(G4int imate, G4String name, G4double a, G4double z,
              G4double dens, G4double radl, G4int nwbf, G4double* ubuf);

void G4gsmixt(G4int imate, G4String name, G4double a[], G4double* z,
              G4double dens, G4int nlmat, G4double* wmat);

void G4gstmed(G4int itmed, G4String name, G4int nmat, G4int isvol,
              G4int ifield, G4double fieldm, G4double tmaxfd,
              G4double stemax, G4double deemax, G4double epsil,
              G4double stmin, G4double* par, G4int npar);

void G4gstpar(G4int itmed, G4String chpar, G4double parval);

void G4gspart(G4int ipart, G4String chnpar, G4int itrtyp, G4double amass,
              G4double charge, G4double tlife, G4double* ubuf,
              G4int nwb);

void G4gsdk(G4int ipart, G4double* bratio, G4int* mode);

void G4gsdet(G4String chset, G4String chdet, G4int nv, G4String* chnmsv,
             G4int* nbitsv, G4int idtyp, G4int nwhi, G4int nwdi);

void G4gsdetv(G4String chset, G4String chdet, G4int idtyp, G4int nwhi,
              G4int nwdi);

void G4gsdeta(G4String chset, G4String chdet, G4String chali,
              G4int nwhi, G4int nwdi);

void G4gsdeth(G4String chset, G4String chdet, G4int nh, G4String* chnamh,
              G4int* nbitsh, G4double* orig, G4double* fact);

void G4gsdetd(G4String chset, G4String chdet, G4int nd, G4String* chnmsd,
              G4int* nbitsd);

void G4gsdetu(G4String chset, G4String chdet, G4int nupar, G4double* upar);

void G4ggclos();

G4LogicalVolume* G4BuildGeom(G4String& inFile);



