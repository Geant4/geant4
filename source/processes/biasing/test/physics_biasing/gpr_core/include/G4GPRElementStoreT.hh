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
// $Id: G4GPRElementStoreT.hh,v 1.6 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRELEMENTSTORET_HH
#define G4GPRELEMENTSTORET_HH

#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRProcessLists.hh"
#include "G4GPRTypeList.hh"
#include "G4GPRSeedManagerT.hh"
#include "G4GPRSingleProcessRelayManagerT.hh"
#include "G4GPRMultiProcessRelayManagerT.hh"
#include "G4GPRSeedT.hh"
#include "G4GPRSingleProcessRelayT.hh"
#include "G4GPRMultiProcessRelayT.hh"
#include "G4GPRMask.hh"
#include "G4GPRMaskManagerT.hh"

template <typename List>
struct G4GPRElementStoreT : G4GPRLinearHierarchyT< G4GPRTypeList_4(G4GPRManagerT< G4GPRSeedT<List> >,
								   G4GPRManagerT< G4GPRSingleProcessRelayT<List> >,
								   G4GPRManagerT< G4GPRMultiProcessRelayT<List> >,
								   G4GPRManagerT< G4GPRMask > ) > {};

/*struct G4GPRElementStoreT : G4GPRLinearHierarchyT< G4GPRTypeList_2(G4GPRManagerT< G4GPRSeedT<List> >,
								   G4GPRManagerT< G4GPRSingleProcessRelayT<List> > ) > {};
*/
#endif
