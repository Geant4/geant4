// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiSigmacPlus.cc,v 1.1 1999-01-07 16:09:59 gunter Exp $
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
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4AntiSigmacPlus.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           AntiSigmacPlus                            ###
// ######################################################################

G4AntiSigmacPlus::G4AntiSigmacPlus(
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
 //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[1];
  // anti_sigmac+ -> anti_lambda_c + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("anti_sigma_c+",1.000,2,"anti_lambda_c+","pi0");

  for (G4int index=0; index <1; index++ ) table->Insert(mode[index]);  
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

G4AntiSigmacPlus G4AntiSigmacPlus::theAntiSigmacPlus(
      "anti_sigma_c+",      2.4535*GeV,       0.0*MeV,  -1.0*eplus, 
		    1,              +1,             0,          
		    2,              -1,             0,             
	     "baryon",               0,            -1,       -4212,
		false,          0.0*ns,          NULL
);

G4AntiSigmacPlus* G4AntiSigmacPlus::AntiSigmacPlusDefinition(){return &theAntiSigmacPlus;}
G4AntiSigmacPlus* G4AntiSigmacPlus::AntiSigmacPlus(){return &theAntiSigmacPlus;}
// initialization for static cut values
G4double   G4AntiSigmacPlus::theAntiSigmacPlusLengthCut = -1.0;
G4double*  G4AntiSigmacPlus::theAntiSigmacPlusKineticEnergyCuts = NULL;

