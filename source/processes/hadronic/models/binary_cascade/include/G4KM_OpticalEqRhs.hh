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
//      GEANT 4 class header file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4KM_OpticalEqRhs.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4KM_OpticalEqRhs_hh
#define G4KM_OpticalEqRhs_hh

#include "globals.hh"
//#include "G4EquationOfMotion.hh"
#include "G4Mag_EqRhs.hh"
#include "G4KM_DummyField.hh" // needed by G4Mag_EqRhs constructor.
#include "G4V3DNucleus.hh"

//class G4KM_OpticalEqRhs : public G4EquationOfMotion
class G4KM_OpticalEqRhs : public G4Mag_EqRhs
{

public:
  G4KM_OpticalEqRhs(G4KM_DummyField *field, G4V3DNucleus * nucleus);
  ~G4KM_OpticalEqRhs();

  virtual void EvaluateRhsGivenB(const G4double y[], const G4double B[3],
				 G4double dydx[]) const;
  virtual void SetChargeMomentumMass(G4ChargeState particleCharge,
				     G4double MomentumXc,
				     G4double MassXc2);
  void SetFactor(G4double mass, G4double opticalParameter);

private:
  G4V3DNucleus * theNucleus;
  G4double theFactor;
  G4double theMass;
};

inline G4KM_OpticalEqRhs::~G4KM_OpticalEqRhs()
{ }


#endif



