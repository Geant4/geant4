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
// $Id: G4GeometryCell.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4GeometryCell
//
// Class description:
//
// This class is usde by scoring and importance sampling.
// It serves to address a "cell". A "cell" is somewhat
// related to a touchable in Geant4. It is identified by a reference
// to a G4VPhysicalVolume and a number (replica number).
// Only simple replicas are supported.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCell_hh
#define G4GeometryCell_hh G4GeometryCell_hh 

#include "globals.hh"

class G4VPhysicalVolume;

class G4GeometryCell
{
public:  // with description

  G4GeometryCell(const G4VPhysicalVolume &aVolume, G4int RepNum);

  G4GeometryCell(const G4GeometryCell &rhs);
    // initialise volume and replica number

  ~G4GeometryCell();
    // simple destruction
  
  const G4VPhysicalVolume &GetPhysicalVolume() const;
    // return the physical volume of the cell

  G4int GetReplicaNumber() const;
    // returns the replica number of the cell

  G4GeometryCell &operator=(const G4GeometryCell &rhs);
  
private:
  const G4VPhysicalVolume *fVPhysicalVolume;
    // pointer to the G4VPhysicalVolume of the "cell" 
    // it is treated as identifyer 

  G4int fRepNum;
    // replica number of the "cell"
};

G4bool operator==(const G4GeometryCell &k1, const G4GeometryCell &k2);
G4bool operator!=(const G4GeometryCell &k1, const G4GeometryCell &k2);

#endif
