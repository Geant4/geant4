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
// $Id: ExN05DetectorConstruction.hh,v 1.4 2002-01-09 17:24:18 ranjard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef ExN05DetectorConstruction_h
#define ExN05DetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4UserLimits;
class ExN05DetectorMessenger;

class ExN05DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ExN05DetectorConstruction();
  ~ExN05DetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();

  
  // methods for UserLimits in crystal
  void      UseUserLimits(G4bool value); 
  G4bool    IsUseUserLimits() { return fUseUserLimits; } 
  G4double  GetMaxTimeInCrystal() const  { return theMaxTimeCutsInCrystal;  }
  G4double  GetMinEkineInCrystal() const { return theMinEkineCutsInCrystal; }
  G4double  GetMinRangeInCrystal() const { return theMinRangeCutsInCrystal; }
  void  SetMaxTimeInCrystal(G4double value);  
  void  SetMinEkineInCrystal(G4double value);  
  void  SetMinRangeInCrystal(G4double value); 

private:
  G4LogicalVolume* theCrystalLog;
  G4LogicalVolume* theTowerLog;

  G4bool           fUseUserLimits;
  G4UserLimits*    theUserLimitsForCrystal; 
  G4double         theMaxTimeCutsInCrystal;
  G4double         theMinEkineCutsInCrystal;
  G4double         theMinRangeCutsInCrystal;

  // messeneger
  ExN05DetectorMessenger* theMessenger;
  
};


#endif


