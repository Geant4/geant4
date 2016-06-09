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
// $Id: G4MassGCellFinder.hh,v 1.3 2006/06/29 21:09:52 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4ParallelGCellFinder
//
// Class description:
//
// Finds a G4GeometryCell in the mass geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassGCellFinder_hh
#define G4MassGCellFinder_hh G4MassGCellFinder_hh

#include "G4Types.hh"
#include "G4VGCellFinder.hh"

class  G4MassGCellFinder : public G4VGCellFinder 
{

public:  // with description

  G4MassGCellFinder();
  virtual  ~G4MassGCellFinder();

  virtual G4GeometryCell 
  GetPreGeometryCell(const G4Step &aStep) const;
  virtual G4GeometryCell 
  GetPostGeometryCell(const G4Step &aStep) const;

};

#endif
