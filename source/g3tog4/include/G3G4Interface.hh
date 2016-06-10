//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G3G4Interface.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
//   Interfaces for G3 equivalent routines
//

#ifndef G3G4INTERFACE_HH
#define G3G4INTERFACE_HH 1

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

void G4gsbool(G4String volName, G4String manyVolName);

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
#endif







