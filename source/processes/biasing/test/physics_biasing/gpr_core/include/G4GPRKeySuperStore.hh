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
// $Id: G4GPRKeySuperStore.hh,v 1.4 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRKEYSUPERSTORE_HH
#define G4GPRKEYSUPERSTORE_HH

#include "G4GPRAssocT.hh"
#include "G4GPRSingletonHierarchyT.hh"
#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRKeyManagerT.hh"
#include "G4GPRProcessLists.hh"

typedef G4GPRLinearHierarchyT< G4GPRTypeList_6(G4GPRKeyManagerT<G4GPRProcessLists::AtRestDoIt>, 
					       G4GPRKeyManagerT<G4GPRProcessLists::AtRestGPIL>,
					       G4GPRKeyManagerT<G4GPRProcessLists::ContinuousDoIt>, 
					       G4GPRKeyManagerT<G4GPRProcessLists::ContinuousGPIL>,
					       G4GPRKeyManagerT<G4GPRProcessLists::DiscreteDoIt>, 
					       G4GPRKeyManagerT<G4GPRProcessLists::DiscreteGPIL>) > G4GPRKeyStore;

typedef G4GPRAssocT<G4GPRPhysicsList*, G4GPRKeyStore> G4GPRPhysicsListAndKey;
typedef G4GPRAssocT<G4ParticleDefinition*, G4GPRPhysicsListAndKey> G4GPRParticleAndKey;

typedef G4GPRSingletonHierarchyT< G4GPRTypeList_1(G4GPRParticleAndKey) > G4GPRKeySuperStore;


#endif
