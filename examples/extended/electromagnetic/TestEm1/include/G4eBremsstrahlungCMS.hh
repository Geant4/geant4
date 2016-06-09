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
// $Id: G4eBremsstrahlungCMS.hh,v 1.3 2004/11/01 09:57:11 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremsstrahlungCMS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 12.09.03
//
// Modifications:
//
// 17-10-03 model variant (V.Ivanchenko)
//
// Class Description:
//
// This class manages the bremsstrahlungCMS for e-/e+
// If energy of emitted gamma is higher than threshold,
// then incoming particle is killed and 2 secondaries are produced 
//

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlungCMS_h
#define G4eBremsstrahlungCMS_h 1

#include "G4eBremsstrahlung.hh"

class G4eBremsstrahlungCMS : public G4eBremsstrahlung
{

public:

  G4eBremsstrahlungCMS(const G4String& name = "eBrem", G4double thresh=DBL_MAX);

 ~G4eBremsstrahlungCMS();

  virtual void SecondariesPostStep( G4VEmModel*,
                              const G4MaterialCutsCouple*,
                              const G4DynamicParticle*,
                                    G4double&,
                                    G4double&);

  virtual void PrintInfoDefinition();

  void SetGammaThreshold(G4double val);

private:

  // hide assignment operator
  G4eBremsstrahlungCMS & operator=(const G4eBremsstrahlungCMS &right);
  G4eBremsstrahlungCMS(const G4eBremsstrahlungCMS&);

  G4double gammaThreshold;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eBremsstrahlungCMS::SetGammaThreshold(G4double val) 
{
  gammaThreshold = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

