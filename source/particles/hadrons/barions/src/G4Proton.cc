// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Proton.cc,v 1.1 1999-01-07 16:10:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  Added particle definitions, H.Kurashige, 19 April 1996
//  Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban, 26/6/96
//  Add ProtonDefinition(), H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4Proton.hh"

// ######################################################################
// ###                           PROTON                               ###
// ######################################################################

G4Proton::G4Proton(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VBarion( aName,mass,width,charge,iSpin,iParity,
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

G4Proton G4Proton::theProton(
             "proton",   0.9382723*GeV,       0.0*MeV,       eplus, 
		    1,              +1,             0,          
		    1,              +1,             0,             
	     "baryon",               0,            +1,        2212,
		 true,            -1.0,          NULL
);

G4Proton* G4Proton::ProtonDefinition(){return &theProton;}
// initialization for static cut values
G4double   G4Proton::theProtonLengthCut = -1.0;
G4double*  G4Proton::theProtonKineticEnergyCuts = NULL;
