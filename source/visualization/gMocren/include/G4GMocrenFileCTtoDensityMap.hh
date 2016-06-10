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
// $Id: G4GMocrenFileCTtoDensityMap.hh 70733 2013-06-05 10:05:57Z gcosmo $
//
//
// Created:  Oct. 12, 2009  Akinori Kimura  
//
// mapping data of  CT values [H.U.] to densities [g/cm3]
//
#ifndef G4GMocrenFile_CTtoDensity_MAP_HH
#define G4GMocrenFile_CTtoDensity_MAP_HH

#include "globals.hh"

class G4GMocrenFileCTtoDensityMap {
public:
  G4GMocrenFileCTtoDensityMap();
  ~G4GMocrenFileCTtoDensityMap();

  // get min. CT value
  G4int GetMinCT() const {return kCTMinMax[0];}
  // get max. CT value
  G4int GetMaxCT() const {return kCTMinMax[1];}
  // get density corresponding to CT value
  G4double GetDensity(G4int & _ct) const;

protected:
  G4int kCTMinMax[2];
  G4double * kDensity;
  G4int kSize;

private:
  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps Coverity happy.
  G4GMocrenFileCTtoDensityMap (const G4GMocrenFileCTtoDensityMap&);
  G4GMocrenFileCTtoDensityMap& operator = (const G4GMocrenFileCTtoDensityMap&);
};
#endif
