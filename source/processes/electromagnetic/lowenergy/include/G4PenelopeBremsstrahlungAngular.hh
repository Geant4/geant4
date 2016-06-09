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
// $Id: G4PenelopeBremsstrahlungAngular.hh,v 1.2 2003/11/07 12:26:08 pandola Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Author: L.Pandola
//
// History:
// -----------
// 04 Feb 2003  L. Pandola       1st implementation
// 07 Nov 2003  L. Pandola       Added method for testing
// Class description:
// Calculation of angular distribution for Penelope Bremsstrahlung
// --------------------------------------------------------------


#ifndef G4PENELOPEBREMSSTRAHLUNGANGULAR_HH
#define G4PENELOPEBREMSSTRAHLUNGANGULAR_HH 1
#include "globals.hh"

class G4PenelopeBremsstrahlungAngular
{ 

private:
  enum{NumberofZPoints=6,
	 NumberofEPoints=6,
	 NumberofKPoints=4,
	 reducedEnergyGrid=21};

public:
 
  G4PenelopeBremsstrahlungAngular(G4int Zmat); 
  ~G4PenelopeBremsstrahlungAngular();
  G4double ExtractCosTheta(G4double PrimaryEnergy,G4double GammaEnergy);
  G4int GetAtomicNumber(); //testing purpose
  
private:

  void InterpolationTableForZ(); //Initialization of tables (part 1)
  void InterpolationForK();  //Initialization of tables (part 2)

  G4double betas[NumberofEPoints]; //betas for interpolation
  //tables for interpolation
  G4double Q1[NumberofEPoints][NumberofKPoints],Q2[NumberofEPoints][NumberofKPoints];
  //expanded tables for interpolation
  G4double Q1E[NumberofEPoints][reducedEnergyGrid],Q2E[NumberofEPoints][reducedEnergyGrid]; 
  //Z of the element
  G4int Zmat;
};


  
#endif
