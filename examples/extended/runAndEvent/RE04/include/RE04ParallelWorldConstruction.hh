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
/// \file runAndEvent/RE04/include/RE04ParallelWorldConstruction.hh
/// \brief Definition of the RE04ParallelWorldConstruction class
//
//
#ifndef RE04ParallelWorldConstruction_h
#define RE04ParallelWorldConstruction_h 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

//
/// A parallel world construction class
///
/// - void Construct()
///     creates a parallel world in the mass world and parameterized volumes
//
class RE04ParallelWorldConstruction : public G4VUserParallelWorld
{
  public:
    RE04ParallelWorldConstruction(G4String& parallelWorldName);
    virtual ~RE04ParallelWorldConstruction();

  public:
    virtual void Construct();

  private:
    G4bool fConstructed;
};

#endif

