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
// $Id: G4ASCIITreeSceneHandler.hh,v 1.15 2006-03-28 17:24:44 allison Exp $
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

#include <set>
#include <vector>
#include <iostream>
#include <fstream>

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
  // of the geometry hierarchy.

  void WriteHeader (std::ostream &);

  const G4VPhysicalVolume* fpLastPV;  // Records last physical volume.
  G4int fPVPCount;              // Counts parameterisations.
  std::ostream* fpOutFile;      // Pointer to output file.
  std::ofstream fOutFile;       // Actual output file (if not G4cout).

  std::set<G4LogicalVolume*> fLVSet;
  typedef std::set<G4LogicalVolume*>::iterator LVSetIterator;
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  std::set<PVPath> fReplicaSet;
  typedef std::set<PVPath>::iterator ReplicaSetIterator;
  typedef std::set<PVPath>::reverse_iterator ReplicaSetReverseIterator;
};

#endif
