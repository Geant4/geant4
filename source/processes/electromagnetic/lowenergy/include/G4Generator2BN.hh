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
// File name:  G4Generator2BN
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta      (luis@lip.pt)
// 
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003                           First implementation acording with new design
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation - 2BN Distribution
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4Generator2BN_h
#define G4Generator2BN_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VBremAngularDistribution.hh"

class G4Generator2BN : public G4VBremAngularDistribution
{

public:

  G4Generator2BN(const G4String& name);

  ~G4Generator2BN();

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
  G4double* Atab;
  G4double* ctab;
  G4double kmin, Ekmin;
  G4double dtheta;
  G4double kcut;

  // hide assignment operator 
  G4Generator2BN & operator=(const  G4Generator2BN &right);
  G4Generator2BN(const  G4Generator2BN&);

};

#endif

