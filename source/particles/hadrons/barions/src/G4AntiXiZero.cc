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
// $Id: G4AntiXiZero.cc,v 1.11 2001-10-24 10:04:59 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  Added particle definitions, H.Kurashige, 14 Feb 1997
// ----------------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

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
 : G4VBaryon( aName,mass,width,charge,iSpin,iParity,
              iConjugation,iIsospin,iIsospin3,gParity,pType,
              lepton,baryon,encoding,stable,lifetime,decaytable )
{
   SetParticleSubType("xi");
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
           "anti_xi0",      1.3149*GeV,  2.27e-12*MeV,         0.0, 
		    1,              +1,             0,          
		    1,              -1,             0,             
	     "baryon",               0,            -1,       -3322,
		false,        0.290*ns,          NULL
);

G4AntiXiZero* G4AntiXiZero::AntiXiZeroDefinition(){return &theAntiXiZero;}
G4AntiXiZero* G4AntiXiZero::AntiXiZero(){return &theAntiXiZero;}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4AntiXiZero::CalcEnergyCuts()
{
  

  // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(LowestEnergy);
}

