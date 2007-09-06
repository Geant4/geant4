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
// $Id: G4GPRTriggerSuperStore.hh,v 1.5 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRTRIGGERSUPERSTORE_HH
#define G4GPRTRIGGERSUPERSTORE_HH

#include "G4GPRPhysicsList.hh"
#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRSingletonHierarchyT.hh"
#include "G4GPRTriggerManagerT.hh"
#include "G4GPRTriggering.hh"

typedef G4GPRLinearHierarchyT< G4GPRTypeList_7(G4GPRTriggerManagerT<G4GPRTriggering::Tracking::StartTracking>, 
					       G4GPRTriggerManagerT<G4GPRTriggering::Tracking::EndTracking>,
					       G4GPRTriggerManagerT<G4GPRTriggering::Stepping::StartStep>, 
					       G4GPRTriggerManagerT<G4GPRTriggering::Stepping::EndStep>,
					       G4GPRTriggerManagerT<G4GPRTriggering::Geometry::StartBoundary>, 
					       G4GPRTriggerManagerT<G4GPRTriggering::Geometry::EndBoundary>,
					       G4GPRTriggerManagerT<G4GPRTriggering::Geometry::NewRegion>) > G4GPRTriggerStore;

typedef G4GPRAssocT<G4GPRPhysicsList*, G4GPRTriggerStore> G4GPRPhysicsListAndTriggerAssoc;
typedef G4GPRAssocT<G4ParticleDefinition*, G4GPRPhysicsListAndTriggerAssoc> G4GPRParticleAndTriggerAssoc;

typedef G4GPRSingletonHierarchyT< G4GPRTypeList_1(G4GPRParticleAndTriggerAssoc) > G4GPRTriggerSuperStore;

#endif
