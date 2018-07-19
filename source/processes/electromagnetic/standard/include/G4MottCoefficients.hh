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
//	G4MottCoefficients.hh
//-------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4MottCoefficients
//
// Author:      Cristina Consolandi
//
// Creation date: 27.05.2012
//
//  Class Description:
//
//       Mott Coulomb Cross section coefficients:
//
//       Reference:
//       M.J. Boschini et al.
//       "Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//       Proc. of the 13th International Conference on Particle Physics and Advanced Technology
//       (13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//
//       Available at: http://arxiv.org/abs/1111.4042v4
//       coeffb of par. 2.1 , eq. (17) were recalculated by M. Tacconi
//       following the same procedur as:
//
//       T. Lijian et al. "Analytic Fitting to the Mott Cross Section of Electrons"
//       Radiat. Phys. Chem. 45 (1995), 235–245.
//
//
// ----------------------------------------------------------------------------------------

//
#ifndef G4MottCoefficients_h
#define G4MottCoefficients_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Pow;

class G4MottCoefficients
{



public:

  	explicit G4MottCoefficients();

  	virtual ~G4MottCoefficients();

	void  SetMottCoeff( G4double targetZ, G4double coeff[5][6] );
        G4double  GetTransitionRandom(G4double targetZ, G4double energy);

private:

 	G4MottCoefficients & operator=(const  G4MottCoefficients &right) = delete;
  	G4MottCoefficients(const  G4MottCoefficients&) = delete;

	G4Pow*          fG4pow;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
