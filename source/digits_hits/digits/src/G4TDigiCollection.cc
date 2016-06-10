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
// $Id: G4TDigiCollection.cc 67992 2013-03-13 10:59:57Z gcosmo $
//

#include "G4TDigiCollection.hh"

G4ThreadLocal G4Allocator<G4DigiCollection> *aDCAllocator_G4MT_TLS_ = 0;

G4DigiCollection::G4DigiCollection() : theCollection((void*)0)
{ if (!aDCAllocator_G4MT_TLS_) aDCAllocator_G4MT_TLS_ = new G4Allocator<G4DigiCollection>  ;;}

G4DigiCollection::G4DigiCollection(G4String detName,G4String colNam)
: G4VDigiCollection(detName,colNam), theCollection((void*)0)
{ if (!aDCAllocator_G4MT_TLS_) aDCAllocator_G4MT_TLS_ = new G4Allocator<G4DigiCollection>  ;;}

G4DigiCollection::~G4DigiCollection()
{ if (!aDCAllocator_G4MT_TLS_) aDCAllocator_G4MT_TLS_ = new G4Allocator<G4DigiCollection>  ;;}

G4int G4DigiCollection::operator==(const G4DigiCollection &right) const
{ if (!aDCAllocator_G4MT_TLS_) aDCAllocator_G4MT_TLS_ = new G4Allocator<G4DigiCollection>  ; return (collectionName==right.collectionName); }

