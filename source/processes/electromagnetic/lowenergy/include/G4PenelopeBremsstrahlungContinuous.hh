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
// $Id: G4PenelopeBremsstrahlungContinuous.hh,v 1.1 2003-03-13 17:23:39 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L.Pandola
//
// History:
// -----------
// 20 Feb 2003  L. Pandola       1st implementation
// Class description:
// Calculation of continuous energy loss for Penelope Bremsstrahlung
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
 
  G4PenelopeBremsstrahlungContinuous(G4int Zmat,G4double tCut, G4double emin, G4double emax); 
  ~G4PenelopeBremsstrahlungContinuous();
  G4double CalculateStopping(G4double PrimaryEnergy);
 
private:

  void PrepareInterpolationTable();
  void LoadFromFile();
 
  //qui ci vanno le cose da tenere a mente fra l'interpolazione 1 e 2
  G4int Zmat;
  G4double tCut;
  G4double MinE,MaxE;
  G4double DLFC; //needed for calculation of extended energy grid
  G4double Energies[NumberofEPoints];
  G4double ReducedCS[NumberofEPoints][NumberofKPoints];
  G4double TotalCS[NumberofEPoints];
  G4double ExtendedLogEnergy[NumberofExtendedEGrid];
  G4double p0[NumberofExtendedEGrid][NumberofKPoints];
  G4double Pbcut[NumberofExtendedEGrid]; //serve?
};


  
#endif
