// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RhoZero.cc,v 1.1 1999-01-07 16:10:20 gunter Exp $
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

#include "G4RhoZero.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         RhoZero                                ###
// ######################################################################


G4RhoZero::G4RhoZero(
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
  // Anti-particle of RhoZero is RhoZero itself  
  SetAntiPDGEncoding(encoding);

  SetPDGStable(false);
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[1];
  // RhoZero -> pi+ + pi-
  mode[0] = new G4PhaseSpaceDecayChannel("rho0",1.000,2,"pi+","pi-");

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

G4RhoZero G4RhoZero::theRhoZero(
	       "rho0",      0.7685*GeV,     150.7*MeV,         0.0, 
		    2,              -1,            -1,          
		    2,               0,            +1,             
	      "meson",               0,             0,         113,
		false,          0.0*ns,          NULL
);

G4RhoZero*    G4RhoZero::RhoZeroDefinition(){return &theRhoZero;}
G4RhoZero*    G4RhoZero::RhoZero(){return &theRhoZero;}
// initialization for static cut values
G4double   G4RhoZero::theRhoZeroLengthCut = -1.0;
G4double*  G4RhoZero::theRhoZeroKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

void G4RhoZero::SetCuts(G4double aCut)
{
  theCutInMaxInteractionLength = aCut;

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  // Create the vector of cuts in energy
  // corresponding to the stopping range cut
  if(theKineticEnergyCuts) delete theKineticEnergyCuts;
  theKineticEnergyCuts = new G4double [materialTable->length()];

  // Build range vector for every material, convert cut into energy-cut,
  // fill theKineticEnergyCuts and delete the range vector
  for (G4int J=0; J<materialTable->length(); J++)
  {
    G4Material* aMaterial = (*materialTable)[J];
    theKineticEnergyCuts[J] = LowestEnergy;
  }
  theRhoZeroLengthCut = theCutInMaxInteractionLength;  
  theRhoZeroKineticEnergyCuts = theKineticEnergyCuts;
}



