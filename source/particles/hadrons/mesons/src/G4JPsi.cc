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
// $Id: G4JPsi.cc,v 1.7 2001-09-19 11:16:54 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      Created                 Hisaya Kurashige, 16 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4JPsi.hh"

#include "G4DecayTable.hh"

// ######################################################################
// ###                          JPsi                            ###
// ######################################################################

G4JPsi::G4JPsi(
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
  SetParticleSubType("J/psi");
  // Anti-particle of J/Psi is J/Psi itself  
  SetAntiPDGEncoding(encoding);


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
G4JPsi G4JPsi::theJPsi(
	      "J/psi",     3.09688*GeV,     0.087*MeV,          0., 
		    2,              -1,            -1,          
		    0,               0,            -1,             
	      "meson",               0,             0,         443,
		false,          0.0*ns,          NULL
);

G4JPsi*  G4JPsi::JPsiDefinition(){return &theJPsi;}
G4JPsi*  G4JPsi::JPsi(){return &theJPsi;}
// initialization for static cut values
G4double   G4JPsi::theJPsiLengthCut = -1.0;
G4double*  G4JPsi::theJPsiKineticEnergyCuts = NULL;

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy

void G4JPsi::SetCuts(G4double aCut)
{
  theCutInMaxInteractionLength = aCut;

  // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(LowestEnergy);

  theJPsiLengthCut = theCutInMaxInteractionLength;  
  theJPsiKineticEnergyCuts = theKineticEnergyCuts;
}





