// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PionZero.cc,v 1.4 2000-02-27 05:57:46 kurasige Exp $
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
//  Operators (+=, *=, ++, -> etc.) correctly used, P. Urban, 26/6/96
//  Add PionZeroDefinition(), H.Kurashige 4 July 1996
// ----------------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4PionZero.hh"


#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          PIONZERO                              ###
// ######################################################################

G4PionZero::G4PionZero(
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
  SetParticleSubType("pi");
  // Anti-particle of Pi0 is pi0 itself  
  SetAntiPDGEncoding(encoding);

  SetPDGStable(false);

  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

  // create a decay channel
  G4VDecayChannel* mode;
  // pi0 -> gamma + gamma
  mode = new G4PhaseSpaceDecayChannel("pi0",0.988,2,"gamma","gamma");
  table->Insert(mode);
  // pi0 -> gamma + e+ + e-
  mode = new G4DalitzDecayChannel("pi0",0.012,"e-","e+");
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

G4PionZero G4PionZero::thePionZero(
		"pi0",   0.1349764*GeV,   7.8e-06*MeV,         0.0, 
		    0,              -1,            +1,          
		    2,               0,            -1,             
	      "meson",               0,             0,         111,
		false,       8.4e-8*ns,          NULL
);

G4PionZero* G4PionZero::PionZeroDefinition(){return &thePionZero;}
// initialization for static cut values
G4double   G4PionZero::thePionZeroLengthCut = -1.0;
G4double*  G4PionZero::thePionZeroKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4PionZero::SetCuts(G4double aCut)
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
  thePionZeroLengthCut = theCutInMaxInteractionLength;  
  thePionZeroKineticEnergyCuts = theKineticEnergyCuts;
}






