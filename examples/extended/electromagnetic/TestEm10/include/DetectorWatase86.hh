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
/// \file electromagnetic/TestEm10/include/DetectorWatase86.hh
/// \brief Definition of the DetectorWatase86 class
//
//
//
// Setup from Y. Watase et al, NIM A248  (1986) 379-388 (fig.7; Li, e-, 2 Gev/c)

#ifndef DetectorWatase86_h
#define DetectorWatase86_h 1

#include "globals.hh"

#include "RadiatorDescription.hh"

class G4VPhysicalVolume;

class DetectorWatase86 
{
  public:
    DetectorWatase86();
    ~DetectorWatase86();

    // methods
    G4VPhysicalVolume* Construct();

    // get methods
    RadiatorDescription* GetRadiatorDescription() const { return fRadiatorDescription; }

  private:
    // data members
    RadiatorDescription* fRadiatorDescription; 
};

#endif
