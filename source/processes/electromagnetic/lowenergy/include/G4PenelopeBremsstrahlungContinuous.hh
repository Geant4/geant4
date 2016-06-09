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
// $Id: G4PenelopeBremsstrahlungContinuous.hh,v 1.2 2003/03/19 10:28:28 pandola Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
