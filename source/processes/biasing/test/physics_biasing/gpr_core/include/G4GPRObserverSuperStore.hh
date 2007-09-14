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
// $Id: G4GPRObserverSuperStore.hh,v 1.1 2007-09-14 16:44:29 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPROBSERVERSUPERSTORE_HH
#define G4GPROBSERVERSUPERSTORE_HH

#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRSingletonHierarchyT.hh"
#include "G4GPRObserverCollectionT.hh"
#include "G4GPRTriggerTypes.hh"
#include "G4GPRObserverT.hh"

typedef G4GPRLinearHierarchyT< G4GPRTypeList_6(G4GPRObserverT<G4GPRTriggerTypes::Initialisation::BuildPhysicsTable>,
					       G4GPRObserverT<G4GPRTriggerTypes::Initialisation::PreparePhysicsTable>,
					       G4GPRObserverT<G4GPRTriggerTypes::Initialisation::RetrievePhysicsTable>,
					       G4GPRObserverT<G4GPRTriggerTypes::Tracking::StartTracking>,
					       G4GPRObserverT<G4GPRTriggerTypes::Tracking::EndTracking>,
					       G4GPRObserverT<G4GPRTriggerTypes::Stepping::StartStep>) > G4GPRObserverStore;

typedef G4GPRAssocT<G4ParticleDefinition*, G4GPRObserverStore> G4GPRParticleAndObserverCollectionStore;

typedef G4GPRSingletonHierarchyT< G4GPRTypeList_1(G4GPRParticleAndObserverCollectionStore) > G4GPRObserverSuperStore;

#endif
