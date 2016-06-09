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
// $Id: G4Penelope01BremsstrahlungAngular.hh,v 1.3 2006-06-29 19:36:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L.Pandola
//
// History:
// -----------
// 04 Feb 2003  L. Pandola       1st implementation
// 07 Nov 2003  L. Pandola       Added method for testing
// 24 May 2011  L. Pandola       Renamed to Penelope01 (obsolete version)
//
// Class description:
// Calculation of angular distribution for Penelope Bremsstrahlung,
// version 2001
// --------------------------------------------------------------


#ifndef G4PENELOPE01BREMSSTRAHLUNGANGULAR_HH
#define G4PENELOPE01BREMSSTRAHLUNGANGULAR_HH 1
#include "globals.hh"

class G4Penelope01BremsstrahlungAngular
{ 

private:
  enum{NumberofZPoints=6,
	 NumberofEPoints=6,
	 NumberofKPoints=4,
	 reducedEnergyGrid=21};

public:
 
  G4Penelope01BremsstrahlungAngular(G4int Zmat); 
  ~G4Penelope01BremsstrahlungAngular();
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
