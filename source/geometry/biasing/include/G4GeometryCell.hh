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
// $Id: G4GeometryCell.hh,v 1.3 2002-09-02 15:22:31 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCell_hh
#define G4GeometryCell_hh G4GeometryCell_hh 

#include "globals.hh"

class G4VPhysicalVolume;

class G4GeometryCell
{
public:
  G4GeometryCell(const G4VPhysicalVolume &aVolume, G4int RepNum);
    // initialise volume and replica number

  ~G4GeometryCell();
    // simple destruction
  
  G4GeometryCell(const G4GeometryCell &rhs);
  G4GeometryCell &operator=(const G4GeometryCell &rhs);

  const G4VPhysicalVolume &GetPhysicalVolume() const
  {return *fVPhysiclaVolume;}

  G4int GetReplicaNumber() const {return fRepNum;}
  
private:
  const G4VPhysicalVolume *fVPhysiclaVolume;
    // pinter to the G4VPhysicalVolume of the "cell" 
    // it is treated as identifyer 
  G4int fRepNum;
    // replica number of the "cell"
};

// -----------------------------------------------------------------------

class G4GeometryCellComp
{

public:  // without description

  G4bool operator() (const G4GeometryCell &k1,
                     const G4GeometryCell &k2) const;
};

// -----------------------------------------------------------------------

G4bool operator==(const G4GeometryCell &k1, const G4GeometryCell &k2);
G4bool operator!=(const G4GeometryCell &k1, const G4GeometryCell &k2);

#endif
