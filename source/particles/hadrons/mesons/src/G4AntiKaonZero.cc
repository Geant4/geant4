// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiKaonZero.cc,v 1.1 1999-01-07 16:10:14 gunter Exp $
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
//                              H.Kurashige   7 Jul 96
// **********************************************************************

#include <fstream.h>
#include <iomanip.h>

#include "G4AntiKaonZero.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                      ANTIKAONZERO                              ###
// ######################################################################

G4AntiKaonZero::G4AntiKaonZero(
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
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[2];
  // kaon0 -> Kaon0L
  mode[0] = new G4PhaseSpaceDecayChannel("anti_kaon0",0.500,1,"kaon0L");
  // kaon0 -> Kaon0S
  mode[1] = new G4PhaseSpaceDecayChannel("anti_kaon0",0.500,1,"kaon0S");

  for (G4int index=0; index <2; index++ ) table->Insert(mode[index]);  
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

G4AntiKaonZero G4AntiKaonZero::theAntiKaonZero(
	 "anti_kaon0",    0.497672*GeV,       0.0*MeV,         0.0, 
		    0,              -1,             0,          
		    1,               1,             0,             
	      "meson",               0,             0,         -311,
		false,             0.0,          NULL
);

G4AntiKaonZero* G4AntiKaonZero::AntiKaonZeroDefinition(){return &theAntiKaonZero;}
// initialization for static cut values
G4double   G4AntiKaonZero::theAntiKaonZeroLengthCut = -1.0;
G4double*  G4AntiKaonZero::theAntiKaonZeroKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy

void G4AntiKaonZero::SetCuts(G4double aCut)
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
  theAntiKaonZeroLengthCut = theCutInMaxInteractionLength;  
  theAntiKaonZeroKineticEnergyCuts = theKineticEnergyCuts;
}

