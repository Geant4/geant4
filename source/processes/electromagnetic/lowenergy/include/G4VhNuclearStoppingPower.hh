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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VhNuclearStoppingPower
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Hadrons/ions nuclear stopping power parameterisation
// Virtual class to provide an interface before nucleare stopping
// power model and G4hLowEnergyIonisation.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4VhNuclearStoppingPower_h
#define G4VhNuclearStoppingPower_h 1

#include "globals.hh"

class G4VhNuclearStoppingPower 
{

public:

  G4VhNuclearStoppingPower();

  virtual ~G4VhNuclearStoppingPower();

  void SetNuclearStoppingFluctuationsOn() {lossFlucFlag = true;}; 

  void SetNuclearStoppingFluctuationsOff() {lossFlucFlag = false;}; 

  virtual G4double NuclearStoppingPower(G4double kineticEnergy,
                                        G4double z1, G4double z2, 
                                        G4double m1, G4double m2) const=0;
protected:

  G4bool lossFlucFlag;
 
private:

  // hide assignment operator 
   G4VhNuclearStoppingPower & operator=(const G4VhNuclearStoppingPower &right);
   G4VhNuclearStoppingPower(const G4VhNuclearStoppingPower&);

};

#endif


