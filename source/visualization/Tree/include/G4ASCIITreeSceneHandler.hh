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
// $Id: G4ASCIITreeSceneHandler.hh,v 1.10 2001-08-24 20:41:26 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy to standard output as
//   ASCII stream.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#ifndef G4ASCIITREESCENEHANDLER_HH
#define G4ASCIITREESCENEHANDLER_HH

#include "G4VTreeSceneHandler.hh"

#include "g4std/set"

class G4VPhysicalVolume;

class G4ASCIITreeSceneHandler: public G4VTreeSceneHandler {
public:
  G4ASCIITreeSceneHandler(G4VGraphicsSystem& system,
			  const G4String& name);
  virtual ~G4ASCIITreeSceneHandler();
  virtual void BeginModeling();
  virtual void EndModeling();
protected:
  virtual void RequestPrimitives(const G4VSolid&);
  // Overrides G4VScenehandler::RequestPrimitives and implements dump
  // of leaves of the geometry heirachy.
  G4std::set<G4LogicalVolume*,G4std::less<G4LogicalVolume*> > fLVSet;
  typedef
  G4std::set<G4LogicalVolume*,G4std::less<G4LogicalVolume*> >::iterator
  LVSetIterator;
  G4std::set<G4VPhysicalVolume*,G4std::less<G4VPhysicalVolume*> > fReplicaSet;
  typedef
  G4std::set<G4VPhysicalVolume*,G4std::less<G4VPhysicalVolume*> >::iterator
  ReplicaSetIterator;
};

#endif
