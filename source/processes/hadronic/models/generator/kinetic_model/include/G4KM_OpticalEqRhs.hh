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
  virtual void SetChargeMomentumMass(G4double particleCharge,
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


// Here by design, but it is unnecessary for nuclear fields
inline void G4KM_OpticalEqRhs::SetChargeMomentumMass(G4double particleCharge,
						     G4double MomentumXc,
						     G4double MassXc2)
{ }

#endif



