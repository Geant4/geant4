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
///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Ecal.hh
// Description: Equipped to describe crystal matrix for different testbeam run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Ecal_h
#define CCalG4Ecal_h 1

#include "CCalEcal.hh"
#include "CCalG4Able.hh"
#include <vector>

typedef G4LogicalVolume* ptrG4Log;

class CCalG4Ecal: public CCalEcal, public CCalG4Able {
public:
  //Backward or Forward type
  enum CMType {module1, module2};

  //Constructor and Destructor
  CCalG4Ecal(const G4String &name);
  virtual ~CCalG4Ecal();

  void setType(CMType ty)    {type = ty;}
  
  //Prefix to all names in the Detector
  static G4String idName;  

protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() ;

private:
  //Methods to construct the different parts of the detector
  G4LogicalVolume* constructGlobal();

private:
  //Private data members
  CMType type;

  //Static logical volumes shared by forward and backward detectors.
  static G4LogicalVolume* crystalmatrixLog;

  // Logical volumes for sensitive detectors
  std::vector<ptrG4Log> sensitiveLogs;

};

#endif
