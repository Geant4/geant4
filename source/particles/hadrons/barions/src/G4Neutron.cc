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
// $Id: G4Neutron.cc,v 1.12 2001-10-24 10:05:00 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//                          H.Kurashige 7 July 1996
//      add neutron life time    Oct 17 2000 
// **********************************************************************
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Neutron.hh"
#include "G4NeutronBetaDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           NEUTRON                              ###
// ######################################################################

G4Neutron::G4Neutron(
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
   SetParticleSubType("nucleon");

   //create Decay Table 
   G4DecayTable*   table = GetDecayTable();
   if (table!=NULL) delete table;
   table = new G4DecayTable();
   
   // create a decay channel
   G4VDecayChannel* mode = new G4NeutronBetaDecayChannel("neutron",1.00);
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

G4Neutron G4Neutron::theNeutron(
            "neutron",  0.93956563*GeV,       0.0*MeV,         0.0, 
		    1,              +1,             0,          
		    1,              -1,             0,             
	     "baryon",               0,            +1,        2112,
		 true,    886.7*second,          NULL
);

G4Neutron* G4Neutron::NeutronDefinition(){return &theNeutron;}
G4Neutron* G4Neutron::Neutron(){return &theNeutron;}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************
//  In this version Input Cut Value is meaning less
//  theKineticEnergyCuts for all materials are set to LowestEnergy
void G4Neutron::CalcEnergyCuts()
{
  

  // Set Energy Cut values to lowest  for all materials
  SetEnergyCutValues(LowestEnergy);

}







