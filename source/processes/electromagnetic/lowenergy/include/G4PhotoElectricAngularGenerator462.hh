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
// File name:  G4PhotoElectricAngularGenerator462
//
// Author:        Andreia Trindade (andreia@lip.pt)
// 
// Creation date: 10 May 2004
//
// Modifications: 
// 10 May 2003       A. Trindade    First implementation acording with new design
//
// Class Description: 
//
// Abstract class for PhotoElectricsstrahlung Angular Generator462 Generation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4PhotoElectricAngularGenerator462_h
#define G4PhotoElectricAngularGenerator462_h 1

#include "G4VPhotoElectricAngularDistribution.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4PhotoElectricAngularGenerator462 : public G4VPhotoElectricAngularDistribution
{

public:

  G4PhotoElectricAngularGenerator462(const G4String& name);

  ~G4PhotoElectricAngularGenerator462();

  G4ThreeVector GetPhotoElectronDirection(G4ThreeVector direction, G4double kineticEnergy);

  void PrintGeneratorInformation() const;

protected:

private:

  // hide assignment operator 
     G4PhotoElectricAngularGenerator462 & operator=(const  G4PhotoElectricAngularGenerator462 &right);
     G4PhotoElectricAngularGenerator462(const  G4PhotoElectricAngularGenerator462&);

};

#endif

