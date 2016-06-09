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
// $Id: G4PenelopeBremsstrahlungContinuous.hh,v 1.3 2006/06/29 19:36:19 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// Author: L.Pandola
//
// History:
// -----------
// 20 Feb 2003  L. Pandola       1st implementation
// 17 Mar 2003  L. Pandola       Added the correction for positrons
// Class description:
// Calculation of continuous energy loss for Penelope Bremsstrahlung
// It is used both for electrons and positrons
// --------------------------------------------------------------


#ifndef G4PENELOPEBREMSSTRAHLUNGCONTINUOUS_HH
#define G4PENELOPEBREMSSTRAHLUNGCONTINUOUS_HH 1
#include "globals.hh"

class G4PenelopeBremsstrahlungContinuous
{ 
private:
  enum{ NumberofEPoints=57,
	  NumberofKPoints=32, 
	  NumberofExtendedEGrid=200};  

public:
 
  G4PenelopeBremsstrahlungContinuous(G4int Zmat,G4double tCut, G4double emin, G4double emax,
				     const G4String partName); 
  ~G4PenelopeBremsstrahlungContinuous();
  G4double CalculateStopping(G4double PrimaryEnergy);
 
private:

  void PrepareInterpolationTable();
  void LoadFromFile();
  G4double PositronCorrection(G4double energy); //correction function for positrons
 
  G4int Zmat;
  G4double tCut;
  G4double MinE,MaxE;
  const G4String partName;
  G4double DLFC; //needed for calculation of extended energy grid
  G4double Energies[NumberofEPoints];
  G4double ReducedCS[NumberofEPoints][NumberofKPoints];
  G4double TotalCS[NumberofEPoints];
  G4double ExtendedLogEnergy[NumberofExtendedEGrid];
  G4double p0[NumberofExtendedEGrid][NumberofKPoints];
  //G4double Pbcut[NumberofExtendedEGrid]; //serve?
};


  
#endif
