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
// $Id: Tst01DetectorConstruction.hh,v 1.2 2001-07-11 10:02:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst01DetectorConstruction_h
#define Tst01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst01DetectorConstruction();
    ~Tst01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
};

#endif

