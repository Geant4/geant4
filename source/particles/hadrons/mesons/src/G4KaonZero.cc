//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4KaonZero.cc,v 1.9 2001-10-24 10:05:31 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//                              H.Kurashige   7 Jul 96
// **********************************************************************

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4KaonZero.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         KAONZERO                              ###
// ######################################################################

G4KaonZero::G4KaonZero(
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
  SetParticleSubType("kaon");
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[2];
  // kaon0 -> Kaon0L
  mode[0] = new G4PhaseSpaceDecayChannel("kaon0",0.500,1,"kaon0L");
  // kaon0 -> Kaon0S
  mode[1] = new G4PhaseSpaceDecayChannel("kaon0",0.500,1,"kaon0S");

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

G4KaonZero G4KaonZero::theKaonZero(
	      "kaon0",    0.497672*GeV,       0.0*MeV,         0.0, 
		    0,              -1,             0,          
		    1,              -1,             0,             
	      "meson",               0,             0,         311,
		false,             0.0,          NULL
);

G4KaonZero* G4KaonZero::KaonZeroDefinition(){return &theKaonZero;}
G4KaonZero* G4KaonZero::KaonZero(){return &theKaonZero;}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4KaonZero::CalcEnergyCuts()
{
  
  
  // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(LowestEnergy);
}

