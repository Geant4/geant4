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
// $Id: XrayFluoNormalization.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03 
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
