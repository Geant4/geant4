// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4.hh,v 1.6 2000-03-07 10:51:40 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I.Hrivnacova, 27 Sep 99

#ifndef G3TOG4_HH
#define G3TOG4_HH 1

#include "globals.hh"

const char gSeparator('_'); 

extern G4int Ipar[1000];
extern G4double Rpar[1000];
extern G4String Spar[1000];

//extern G4double min(G4double, G4double);

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
