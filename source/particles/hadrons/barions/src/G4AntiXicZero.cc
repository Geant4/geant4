// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiXicZero.cc,v 1.3 1999-10-03 09:13:18 kurasige Exp $
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

#include "G4AntiXicZero.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           AntiXicZero                          ###
// ######################################################################

G4AntiXicZero::G4AntiXicZero(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VBaryon( aName,mass,width,charge,iSpin,iParity,
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

G4AntiXicZero G4AntiXicZero::theAntiXicZero(
         "anti_xi_c0",      2.4703*GeV,    6.7e-9*MeV,         0.0, 
		    1,              +1,             0,          
		    1,              +1,             0,             
	     "baryon",               0,            -1,       -4132,
		false,     0.098e-3*ns,          NULL
);

G4AntiXicZero* G4AntiXicZero::AntiXicZeroDefinition(){return &theAntiXicZero;}
G4AntiXicZero* G4AntiXicZero::AntiXicZero(){return &theAntiXicZero;}
// initialization for static cut values
G4double   G4AntiXicZero::theAntiXicZeroLengthCut = -1.0;
G4double*  G4AntiXicZero::theAntiXicZeroKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4AntiXicZero::SetCuts(G4double aCut)
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
  theAntiXicZeroLengthCut = theCutInMaxInteractionLength;  
  theAntiXicZeroKineticEnergyCuts = theKineticEnergyCuts;
  // Rebuild the physics tables for every process for this particle type
  
}

