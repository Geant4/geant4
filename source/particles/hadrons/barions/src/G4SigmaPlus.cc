// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SigmaPlus.cc,v 1.1 1999-01-07 16:10:02 gunter Exp $
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
//  Added particle definitions, H.Kurashige, 14 Feb 1997
// ----------------------------------------------------------------------

#include <fstream.h>
#include <iomanip.h>

#include "G4SigmaPlus.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           SigmaPlus                            ###
// ######################################################################

G4SigmaPlus::G4SigmaPlus(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable )
 : G4VBarion( aName,mass,width,charge,iSpin,iParity,
              iConjugation,iIsospin,iIsospin3,gParity,pType,
              lepton,baryon,encoding,stable,lifetime,decaytable )
{
  //create Decay Table 
  G4DecayTable*   table = GetDecayTable();
  if (table!=NULL) delete table;
  table = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[2];
  // sigma+ -> proton + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("sigma+",0.516,2,"proton","pi0");
  // sigma+ -> neutron + pi+
  mode[1] = new G4PhaseSpaceDecayChannel("sigma+",0.483,2,"neutron","pi+");
 
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

G4SigmaPlus G4SigmaPlus::theSigmaPlus(
             "sigma+",     1.18937*GeV,       0.0*MeV,       eplus, 
		    1,              +1,             0,          
		    2,              +2,             0,             
	     "baryon",               0,            +1,        3222,
		false,       0.0799*ns,          NULL
);

G4SigmaPlus* G4SigmaPlus::SigmaPlusDefinition(){return &theSigmaPlus;}
// initialization for static cut values
G4double   G4SigmaPlus::theSigmaPlusLengthCut = -1.0;
G4double*  G4SigmaPlus::theSigmaPlusKineticEnergyCuts = NULL;
