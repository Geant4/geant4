// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4XicPlus.cc,v 1.2 1999-06-09 16:08:44 kurasige Exp $
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
//  Added particle definitions, H.Kurashige, 14 June 1997
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
//  Modified PDG encoding           H.Kurashige 24 Sep. 98
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4XicPlus.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           XicPlus                              ###
// ######################################################################

G4XicPlus::G4XicPlus(
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

G4XicPlus G4XicPlus::theXicPlus(
              "xi_c+",      2.4656*GeV,   1.86e-9*MeV,   +1.*eplus, 
		    1,              +1,             0,          
		    1,              +1,             0,             
	     "baryon",               0,            +1,        4232,
		false,       0.35e-3*ns,          NULL
);

G4XicPlus* G4XicPlus::XicPlusDefinition(){return &theXicPlus;}
G4XicPlus* G4XicPlus::XicPlus(){return &theXicPlus;}
// initialization for static cut values
G4double   G4XicPlus::theXicPlusLengthCut = -1.0;
G4double*  G4XicPlus::theXicPlusKineticEnergyCuts = NULL;
