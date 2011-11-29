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
// File name:  G4RDGenerator2BN
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta      (luis@lip.pt)
// 
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003                           First implementation acording with new design
// 07 Nov 2003                           Lockup tables for fast initialization
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - 2BN Distribution
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4RDGenerator2BN_h
#define G4RDGenerator2BN_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4RDVBremAngularDistribution.hh"

class G4RDGenerator2BN : public G4RDVBremAngularDistribution
{

public:

  G4RDGenerator2BN(const G4String& name);

  ~G4RDGenerator2BN();

  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  void PrintGeneratorInformation() const;

public:

  void SetInterpolationThetaIncrement(G4double increment) {dtheta = increment;};
  G4double GetInterpolationThetaIncrement() {return dtheta;};

  void SetGammaCutValue(G4double cutValue) {kcut = cutValue;};
  G4double GetGammaCutValue() {return kcut;};

  void ConstructMajorantSurface();

protected:

  G4double CalculateFkt(G4double k, G4double theta, G4double A, G4double c) const;
  G4double Calculatedsdkdt(G4double kout, G4double theta, G4double Eel) const;
  G4double Generate2BN(G4double Ek, G4double k) const;


private:

  G4double b;
  G4int index_min, index_max;
  G4double kmin, Ekmin;
  G4double dtheta;
  G4double kcut;
  static G4double Atab[320];
  static G4double ctab[320];

  // hide assignment operator 
  G4RDGenerator2BN & operator=(const  G4RDGenerator2BN &right);
  G4RDGenerator2BN(const  G4RDGenerator2BN&);

};

#endif

