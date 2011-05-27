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
// $Id: G4PenelopeBremsstrahlungAngular.hh,v 1.1 2010-12-20 14:11:13 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L.Pandola
//
// History:
// -----------
// 23 Nov 2010  L. Pandola       1st implementation
// 24 May 2011  L. Pandola       Renamed (make default Penelope)
//
// Class description:
// Calculation of angular distribution for Penelope Bremsstrahlung
// version 2008
// --------------------------------------------------------------


#ifndef G4PENELOPEBREMSSTRAHLUNGANGULAR_HH
#define G4PENELOPEBREMSSTRAHLUNGANGULAR_HH 1
#include "globals.hh"
#include <map>

class G4PhysicsTable;

class G4PenelopeBremsstrahlungAngular
{ 

public:
  G4PenelopeBremsstrahlungAngular(); 
  ~G4PenelopeBremsstrahlungAngular();
  G4double SampleCosTheta(G4double Zeq,
			  G4double primaryEnergy,G4double gammaEnergy);
  void ClearTables();
  
private:
  void PrepareInterpolationTables(G4double Zeq);
  
  //Tables containing the Lorentz sampling coefficients 
  //The key is the effective Z of the material
  std::map<G4double,G4PhysicsTable*> *theLorentzTables1;
  std::map<G4double,G4PhysicsTable*> *theLorentzTables2;

  void ReadDataFile();
  G4bool dataRead;
  
  static const G4int NumberofZPoints=6;
  static const G4int NumberofEPoints=6;
  static const G4int NumberofKPoints=4;

  G4double QQ1[NumberofZPoints][NumberofEPoints][NumberofKPoints];
  G4double QQ2[NumberofZPoints][NumberofEPoints][NumberofKPoints];

};


  
#endif
