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


