// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RhoMinus.cc,v 1.1 1999-01-07 16:10:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation   Hisaya Kurashige,  8 June 1998
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ------------------------------------------------------------


#include <fstream.h>
#include <iomanip.h>

#include "G4RhoMinus.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         RhoMinus                                ###
// ######################################################################


G4RhoMinus::G4RhoMinus(
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
  SetPDGStable(false);
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[1];
  // RhoMinus -> pi+ + pi-
  mode[0] = new G4PhaseSpaceDecayChannel("rho+",1.000,2,"pi-","pi0");

  table->Insert(mode[0]);  

  delete [] mode;

  SetDecayTable(table);
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

G4RhoMinus G4RhoMinus::theRhoMinus(
	       "rho0",      0.7685*GeV,     150.7*MeV,    -1.0*eplus, 
		    2,              -1,             0,          
		    2,              -2,            +1,             
	      "meson",               0,             0,        -213,
		false,          0.0*ns,          NULL
);

G4RhoMinus*    G4RhoMinus::RhoMinusDefinition(){return &theRhoMinus;}
G4RhoMinus*    G4RhoMinus::RhoMinus(){return &theRhoMinus;}
// initialization for static cut values
G4double   G4RhoMinus::theRhoMinusLengthCut = -1.0;
G4double*  G4RhoMinus::theRhoMinusKineticEnergyCuts = NULL;

