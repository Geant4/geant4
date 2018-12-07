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
/// \file electromagnetic/TestEm10/include/DetectorBari05.hh
/// \brief Definition of the DetectorBari05 class
//
//
//
// Setup for Bari INFN XTR test beam (~2004) at CERN. With He beam-pipe
// M. Brigida et al, NIM A550 (2005) 157-168
// Runs by : TestEm10 bari05.mac

#ifndef DetectorBari05_h
#define DetectorBari05_h 1

#include "globals.hh"

#include "RadiatorDescription.hh"

class G4VPhysicalVolume;

class DetectorBari05 
{
  public:
    DetectorBari05();
    ~DetectorBari05();

    // methods
    G4VPhysicalVolume* Construct();

    // get methods
    RadiatorDescription* GetRadiatorDescription() const { return fRadiatorDescription; }

  private:
    // data members
    RadiatorDescription* fRadiatorDescription; 
};

#endif
