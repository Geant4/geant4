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
// $Id: G4AntiNeutrinoTau.cc,v 1.8 2001-10-16 08:16:17 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo 
//      by             H.Kurashige,7  July 1996
// **********************************************************************

#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4AntiNeutrinoTau.hh"

// ######################################################################
// ###                   ANTI TAU NEUTRINO                            ###
// ######################################################################

G4AntiNeutrinoTau::G4AntiNeutrinoTau(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VLepton( aName,mass,width,charge,iSpin,iParity,
              iConjugation,iIsospin,iIsospin3,gParity,pType,
              lepton,baryon,encoding,stable,lifetime,decaytable )
{
  SetParticleSubType("tau");
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
//
G4AntiNeutrinoTau G4AntiNeutrinoTau::theAntiNeutrinoTau(
	"anti_nu_tau",          0.0*MeV,       0.0*MeV,         0.0, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",              -1,             0,          -16,
		 true,             0.0,          NULL
);

G4AntiNeutrinoTau* G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()
{
  return &theAntiNeutrinoTau;
}

G4AntiNeutrinoTau* G4AntiNeutrinoTau::AntiNeutrinoTau()
{
  return &theAntiNeutrinoTau;
}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

void G4AntiNeutrinoTau::SetCuts(G4double aCut)
{
  SetCutInMaxInteractionLength( aCut );

   // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(0.0*keV); 

}
