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
// G4RToEConvForPositron
//
// Class description:
//
// This class is a Range to Energy Converter for positron.

// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4RToEConvForPositron_hh
#define G4RToEConvForPositron_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4VRangeToEnergyConverter.hh"

class G4RToEConvForPositron : public G4VRangeToEnergyConverter
{
  public:

    G4RToEConvForPositron();
      // Constructor

  virtual ~G4RToEConvForPositron();
      // Destructor

  protected:

    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy);

    G4double Mass = 0.0;
    G4double Z = -1.0;  
    G4double taul = 0.0;
    G4double ionpot = 0.0;
    G4double ionpotlog = -1.0e-10;
    G4double bremfactor = 0.1;
};


#endif
