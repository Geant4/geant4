// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PionPlus.cc,v 1.1 1999-01-07 16:10:19 gunter Exp $
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
//  Add PionPlusDefinition(), H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4PionPlus.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          PIONPLUS                              ###
// ######################################################################

G4PionPlus::G4PionPlus(
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

  // create a decay channel
  // pi+ -> mu+ + nu_mu
  G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel("pi+",1.00,2,"mu+","nu_mu");
  table->Insert(mode);

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

// In this version, charged pions are set to stable
G4PionPlus G4PionPlus::thePionPlus(
		 "pi+",  0.1395700*GeV,       0.0*MeV,    +1.*eplus, 
		    0,              -1,             0,          
		    2,              +2,            -1,             
	      "meson",               0,             0,         211,
		 true,       26.030*ns,          NULL
);

G4PionPlus*  G4PionPlus::PionPlusDefinition(){return &thePionPlus;}
// initialization for static cut values
G4double   G4PionPlus::thePionPlusLengthCut = -1.0;
G4double*  G4PionPlus::thePionPlusKineticEnergyCuts = NULL;





