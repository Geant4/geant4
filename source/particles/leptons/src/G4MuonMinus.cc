// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuonMinus.cc,v 1.2 1999-06-09 15:37:32 kurasige Exp $
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
//  Added SetCuts implementation, L.Urban, 12 June 1996
//  Operators (+=, *=, ++, -> etc.) correctly used, P. Urban, 26/6/96
//  Add MuonMinusDefinition(), H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4MuonMinus.hh"

#include "G4MuonDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          MUONMINUS                             ###
// ######################################################################

G4MuonMinus::G4MuonMinus(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VLepton( aName,mass,width,charge,iSpin,iParity,
              iConjugation,iIsospin,iIsospin3,gParity,pType,
              lepton,baryon,encoding,stable,lifetime,decaytable )
{
  SetPDGStable(false);
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

  // create a decay channel
  G4VDecayChannel* mode = new G4MuonDecayChannel("mu-",1.00);
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

G4MuonMinus G4MuonMinus::theMuonMinus(
		 "mu-",  0.1056584*GeV, 2.99591e-16*MeV,  -1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",               1,             0,          13,
		 true,      2197.03*ns,          NULL
);

G4MuonMinus* G4MuonMinus::MuonMinusDefinition() {return &theMuonMinus;}
G4double   G4MuonMinus::theMuonMinusLengthCut = -1.0;
G4double*  G4MuonMinus::theMuonMinusKineticEnergyCuts = NULL;
