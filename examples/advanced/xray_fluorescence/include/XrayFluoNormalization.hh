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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------
#ifndef XrayFluoNormalization_h
#define XrayFluoNormalization_h 1
#include "globals.hh"

class XrayFluoDataSet;
class XrayFluoNormalization
{
public:

  //constructor
  XrayFluoNormalization();

  //destructor
  ~XrayFluoNormalization();

  //this method returns a data set equivalent to the one in the file whose 
  //name must be passed as the last argument normalized to the value returned 
  //by Integrate
  //the first and second arguments identifies the energy value in which 
  //Integrate() integrates, the third is the number of buns used in the
  //integration
  const XrayFluoDataSet* Normalize(G4double, G4double, G4int,G4String);

private:
  //this method integrates the function achieved interpolating
  //berween the points of the data file and returns the value of the integral 
  G4double Integrate(G4double, G4double, G4int, XrayFluoDataSet*);

};
 
#endif
