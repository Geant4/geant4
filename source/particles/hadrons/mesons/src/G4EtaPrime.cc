// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EtaPrime.cc,v 1.3 1999-12-15 14:51:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, 8 June 1998 Hisaya Kurashige
// ****************************************************************
// ----------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4EtaPrime.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         EtaPrime                                    ###
// ######################################################################


G4EtaPrime::G4EtaPrime(
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
  // Anti-particle of EtaPrime is EtaPrime itself  
  SetAntiPDGEncoding(encoding);

  SetPDGStable(false);

  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[3];
  // EtaPrime -> eta + pi+ + pi-
  mode[0] = new G4PhaseSpaceDecayChannel("eta_prime",0.437,3,"eta","pi+","pi-");
  // EtaPrime -> eta + pi0 + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("eta_prime",0.208,3,"eta","pi0","pi0");
  // EtaPrime -> rho0 + gamma
  mode[2] = new G4PhaseSpaceDecayChannel("eta_prime",0.302,2,"rho0","gamma");

  for (G4int index=0; index <3; index++ ) table->Insert(mode[index]);  
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

G4EtaPrime G4EtaPrime::theEtaPrime(
	  "eta_prime",     0.95777*GeV,     0.203*MeV,         0.0, 
		    0,              -1,            +1,          
		    0,               0,            +1,             
	      "meson",               0,             0,         331,
		false,          0.0*ns,          NULL
);

G4EtaPrime*    G4EtaPrime::EtaPrimeDefinition(){return &theEtaPrime;}
// initialization for static cut values
G4double   G4EtaPrime::theEtaPrimeLengthCut = -1.0;
G4double*  G4EtaPrime::theEtaPrimeKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

void G4EtaPrime::SetCuts(G4double aCut)
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
  theEtaPrimeLengthCut = theCutInMaxInteractionLength;  
  theEtaPrimeKineticEnergyCuts = theKineticEnergyCuts;
}



