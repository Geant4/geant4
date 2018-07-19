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
// $Id: G4MuonicAtom.cc 98732 2016-08-09 10:50:57Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: 
//      July 2016, K.Lynch - first implementation
//      June 2017, K.L.Genser - added baseion, lifetimes and access functions

#include "G4MuonicAtom.hh"

G4MuonicAtom::G4MuonicAtom(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifeTime,
       G4DecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType,      G4Ions const*       baseion,
       G4int               anti_encoding,
       G4double            excitation, 
       G4int               isomer,
       G4double            DIOLifeTime,
       G4double            NCLifeTime
   ) : G4Ions( aName,mass,width,charge,iSpin,iParity,
               iConjugation,iIsospin,iIsospin3,gParity,pType,
               lepton,baryon,encoding,stable,lifeTime,decaytable,
               shortlived,subType,anti_encoding,excitation,isomer),
       baseIon(baseion),
       fDIOLifeTime(DIOLifeTime),
       fNCLifeTime(NCLifeTime)
{
  // member is private in base, so we need to go through the public
  // interface
  SetFloatLevelBase(G4FloatLevelBase::no_Float);
  // from G4ParticleDefinition
  isGeneralIon = false;
  isMuonicAtom = true;
}

G4MuonicAtom::~G4MuonicAtom(){}
