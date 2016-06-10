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
// $Id: G3toG4.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// modified by I.Hrivnacova, 27 Sep 99

#ifndef G3TOG4_HH
#define G3TOG4_HH 1

#include "G3toG4Defs.hh"
#include "globals.hh"

extern G3G4DLL_API char gSeparator;
extern G3G4DLL_API G4int Ipar[1000];
extern G3G4DLL_API G4double Rpar[1000];
extern G3G4DLL_API G4String Spar[1000];

void G3fillParams(G4String *tokens, const char *ptypes);
// G4bool G3NegVolPars(G4double pars[], G4int* np, G4String vol, G4String moth,
//               char* routine);

#include "G3G4Interface.hh"

// Parameter types for Geant routines
//   s=string  i=integer r=real   capitalized=array
//   In case of arrays, the last integer before the array is the 
//   number of elements.
#define PTgsvolu "ssiiR"
#define PTgspos  "sisrrris"
#define PTgsposp "sisrrrisiR"
#define PTgsatt  "ssi"
#define PTgsrotm "irrrrrr"
#define PTgsdvn  "ssii"
#define PTgsdvt  "ssriii"
#define PTgsdvx  "ssiirrii"
#define PTgsdvn2 "ssiiri"
#define PTgsdvt2 "ssririi"
#define PTgsmate "isrrrriR"
//#define PTgsmixt "isriRRR"
#define PTgsmixt "isriQ"
#define PTgstmed "isiiirrrrrriR"
#define PTgstpar "isr"
#define PTgspart "isirrriR"
#define PTgsdk   "iiRI"
#define PTgsdet  "ssiSIiii"
#define PTgsdetv "ssiii"
#define PTgsdeta "sssii"
#define PTgsdeth "ssiSIRR"
#define PTgsdetd "ssiSI"
#define PTgsdetu "ssiR"
#endif
