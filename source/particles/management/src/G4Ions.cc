//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Ions.cc 98732 2016-08-09 10:50:57Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************

#include <fstream>
#include <iomanip>

#include "G4Ions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

 // ######################################################################
// ###                           Ions                                ###
// ######################################################################

G4Ions::G4Ions(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable , G4bool              shortlived,
       const G4String&     subType,
       G4int               anti_encoding,
       G4double            excitation,
       G4int               isomer
)
  : G4ParticleDefinition( aName,mass,width,charge,iSpin,iParity,
           iConjugation,iIsospin,iIsospin3,gParity,pType,
           lepton,baryon,encoding,stable,lifetime,decaytable,
			  shortlived, subType, anti_encoding),
    theExcitationEnergy(excitation),
    theIsomerLevel(isomer),
    floatLevelBase(G4FloatLevelBase::no_Float)
{
   if ((aName == "proton") || (aName == "neutron")) { 
     isGeneralIon = false ;
   } else if (   (aName == "GenericIon") || (aName == "alpha") 
       || (aName == "He3") || (aName == "deuteron")|| (aName == "triton")) {
     isGeneralIon = false ;
   } else if ( (aName == "anti_He3") || (aName == "anti_deuteron")
	|| (aName == "anti_triton") || (aName == "anti_alpha") ) {
     isGeneralIon = false ;
   } else if ( (aName == "iron") || (aName == "oxygen") ||  (aName == "nitrogen")
        || (aName == "carbon") || (aName == "helium") || (aName == "alpha+")
        || (aName == "hydrogen") || (aName == "Ps-1s") || (aName == "Ps-2s")) {
     isGeneralIon = false ;
   } else {
     isGeneralIon = true;
   }

  // isomer level isset to 9 
  // if isomer level is set to 0 for excited state
  if ((theExcitationEnergy > 0.0) && (isomer==0)) isomer =9; 

   if (GetAtomicNumber() == 0  ) {
     // AtomicNumber/Mass is positve even for anti_nulceus
     SetAtomicNumber( std::abs(G4int(GetPDGCharge()/eplus)) );
     SetAtomicMass( std::abs(GetBaryonNumber()) );
   }
}


G4Ions::~G4Ions()
{
}





