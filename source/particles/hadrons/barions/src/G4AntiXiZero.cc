// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiXiZero.cc,v 1.1 1999-01-07 16:10:00 gunter Exp $
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
//  Added particle definitions, H.Kurashige, 14 Feb 1997
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4AntiXiZero.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           AntiXiZero                           ###
// ######################################################################

G4AntiXiZero::G4AntiXiZero(
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
  // anti_xi0 -> anti_lambda + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("anti_xi0",1.000,2,"anti_lambda","pi0");
 
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

G4AntiXiZero G4AntiXiZero::theAntiXiZero(
           "anti_xi0",      1.3149*GeV,       0.0*MeV,         0.0, 
		    1,              +1,             0,          
		    1,              -1,             0,             
	     "baryon",               0,            -1,       -3322,
		false,        0.290*ns,          NULL
);

G4AntiXiZero* G4AntiXiZero::AntiXiZeroDefinition(){return &theAntiXiZero;}
// initialization for static cut values
G4double   G4AntiXiZero::theAntiXiZeroLengthCut = -1.0;
G4double*  G4AntiXiZero::theAntiXiZeroKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4AntiXiZero::SetCuts(G4double aCut)
{
  theCutInMaxInteractionLength = aCut;

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  // Create the vector of cuts in energy
  // corresponding to the stopping range cut
  if(theKineticEnergyCuts) delete [] theKineticEnergyCuts;
  theKineticEnergyCuts = new G4double [materialTable->length()];

  // Build range vector for every material, convert cut into energy-cut,
  // fill theKineticEnergyCuts and delete the range vector
  for (G4int J=0; J<materialTable->length(); J++)
  {
    G4Material* aMaterial = (*materialTable)[J];
    theKineticEnergyCuts[J] = LowestEnergy;
  }
  theAntiXiZeroLengthCut = theCutInMaxInteractionLength;  
  theAntiXiZeroKineticEnergyCuts = theKineticEnergyCuts;
  // Rebuild the physics tables for every process for this particle type
  
}

