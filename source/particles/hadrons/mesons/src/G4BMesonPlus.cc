// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BMesonPlus.cc,v 1.1 1999-01-07 16:10:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//
//      Created                 Hisaya Kurashige, 16 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4BMesonPlus.hh"

#include "G4DecayTable.hh"

// ######################################################################
// ###                          BMesonPLUS                            ###
// ######################################################################

G4BMesonPlus::G4BMesonPlus(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VMeson( aName,mass,width,charge,iSpin,iParity,
             iConjugation,iIsospin,iIsospin3,gParity,pType,
             lepton,baryon,encoding,stable,lifetime,decaytable )
{
}

// ......................................................................
// ...                 static member definitions                      ...
// ......................................................................
//     
//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 

// In this version, charged pions are set to stable
G4BMesonPlus G4BMesonPlus::theBMesonPlus(
	         "B+",      5.2789*GeV,       0.0*MeV,    +1.*eplus, 
		    0,              -1,             0,          
		    1,              +1,             0,             
	      "meson",               0,             0,         521,
		false,      1.62e-3*ns,          NULL
);

G4BMesonPlus*  G4BMesonPlus::BMesonPlusDefinition(){return &theBMesonPlus;}
G4BMesonPlus*  G4BMesonPlus::BMesonPlus(){return &theBMesonPlus;}
// initialization for static cut values
G4double   G4BMesonPlus::theBMesonPlusLengthCut = -1.0;
G4double*  G4BMesonPlus::theBMesonPlusKineticEnergyCuts = NULL;





