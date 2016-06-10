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
// $Id: G4ASCIITreeSceneHandler.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
#include <sstream>

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

  std::ostream* fpOutFile;      // Pointer to output file.
  std::ofstream fOutFile;       // Actual output file (if not G4cout).
  std::ostringstream fRestOfLine;  // For accumulating rest of line
				   // pending arrival of new
				   // G4PVPlacement.
  const G4VPhysicalVolume* fpLastPV;
  G4String fLastPVName;
  G4int    fLastCopyNo;
  G4int    fLastNonSequentialCopyNo;

  std::set<G4LogicalVolume*> fLVSet;
  typedef std::set<G4LogicalVolume*>::iterator LVSetIterator;
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  std::set<PVPath> fReplicaSet;
  typedef std::set<PVPath>::iterator ReplicaSetIterator;
  typedef std::set<PVPath>::reverse_iterator ReplicaSetReverseIterator;
};

#endif
