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
// $Id: G4AntiBsMesonZero.cc,v 1.8 2001-10-15 10:08:22 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      Created                 Hisaya Kurashige, 16 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4AntiBsMesonZero.hh"

#include "G4DecayTable.hh"

// ######################################################################
// ###                      AntiBsMesonZero                            ###
// ######################################################################

G4AntiBsMesonZero::G4AntiBsMesonZero(
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
   SetParticleSubType("Bs");
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

// In this version, charged pions are set to stable
G4AntiBsMesonZero G4AntiBsMesonZero::theAntiBsMesonZero(
	   "anti_Bs0",      5.3692*GeV,   4.27e-10*MeV,          0., 
		    0,              -1,             0,          
		    0,               0,             0,             
	      "meson",               0,             0,        -531,
		false,      1.61e-3*ns,          NULL
);

G4AntiBsMesonZero*  G4AntiBsMesonZero::AntiBsMesonZeroDefinition(){return &theAntiBsMesonZero;}
G4AntiBsMesonZero*  G4AntiBsMesonZero::AntiBsMesonZero(){return &theAntiBsMesonZero;}

void G4AntiBsMesonZero::SetCuts(G4double aCut)
{
  SetCutInMaxInteractionLength( aCut );

  // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(LowestEnergy);

}






