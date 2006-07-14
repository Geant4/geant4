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
// $Id: A01DetectorConstruction.hh,v 1.1 2006-07-14 14:42:46 asaim Exp $
// --------------------------------------------------------------
//

#ifndef A01DetectorConstruction_h
#define A01DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4Material;
class G4VisAttributes;

class A01DetectorConstruction : public G4VUserDetectorConstruction 
{
public:
    A01DetectorConstruction();
    virtual ~A01DetectorConstruction();

public:
    virtual G4VPhysicalVolume* Construct();

private:
    void ConstructMaterials();

private:
    G4Material* air;
    G4Material* water;

    G4VisAttributes* worldVisAtt;
    G4VisAttributes* phantomVisAtt;

};

#endif

