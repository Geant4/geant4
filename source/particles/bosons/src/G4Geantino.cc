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
// $Id: G4Geantino.cc,v 1.11 2002-12-16 11:15:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  Added particle definitions, H.Kurashige, 19 April 1996
//  Operators (+=, *=, ++, -> etc.) correctly used, P. Urban, 26/6/96
//  Add G4Geantino::GeantinoDefinition() H.Kurashige, 4 July 1996
// ----------------------------------------------------------------------

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Geantino.hh"



// ######################################################################
// ###                          GEANTINO                              ###
// ######################################################################

G4Geantino::G4Geantino(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VBoson( aName,mass,width,charge,iSpin,iParity,
             iConjugation,iIsospin,iIsospin3,gParity,pType,
             lepton,baryon,encoding,stable,lifetime,decaytable )
{
   SetParticleSubType("geantino");
   // Anti-particle of geantino is geantino itself
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

G4Geantino G4Geantino::theGeantino(
	  "geantino",          0.0*MeV,       0.0*MeV,         0.0, 
		    0,               0,             0,          
		    0,               0,             0,             
	   "geantino",               0,             0,           0,
		 true,             0.0,          NULL
);

G4Geantino* G4Geantino::GeantinoDefinition() {return &theGeantino;}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

G4Geantino* G4Geantino::Geantino()
{
  return &theGeantino; 
}
