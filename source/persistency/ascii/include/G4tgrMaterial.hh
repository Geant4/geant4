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
// $Id: G4tgrMaterial.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrMaterial
//
// Class description:
//
// Class to represent a transient material.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------
#ifndef G4tgrMaterial_h
#define G4tgrMaterial_h

#include "globals.hh"
#include "G4Material.hh"

class G4tgrMaterial
{
  public:  // with description

    G4tgrMaterial();
    virtual ~G4tgrMaterial();

    // Get methods

    const G4String& GetName() const { return theName; }  

    G4double GetDensity() const { return theDensity; }
      // Density in g/cm3

    G4int GetNumberOfComponents() const { return theNoComponents; } 
    const G4String& GetType() const { return theMateType; }  

    virtual G4double GetA() const = 0;
    virtual G4double GetZ() const = 0;
    virtual const G4String& GetComponent(G4int i) const = 0;
    virtual G4double GetFraction(G4int i) = 0;

    G4double GetIonisationMeanExcitationEnergy() const
      { return theIonisationMeanExcitationEnergy; }

    void SetIonisationMeanExcitationEnergy( G4double mee )
      { theIonisationMeanExcitationEnergy = mee; }

    G4State GetState() const { return theState; }
    void SetState( G4String val );

    G4double GetTemperature() const { return theTemperature; }
    void SetTemperature( G4double val ) { theTemperature = val; }

    G4double GetPressure() const { return thePressure; }
    void SetPressure( G4double val ) { thePressure = val; }

  protected:

    G4String  theName;
    G4double  theDensity;
    G4int theNoComponents;
    G4String  theMateType;
    G4double theIonisationMeanExcitationEnergy;
    G4State theState;
    G4double theTemperature;
    G4double thePressure;
};

#endif
