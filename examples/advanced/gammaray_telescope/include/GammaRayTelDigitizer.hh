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
// $Id: GammaRayTelDigitizer.hh,v 1.2 2001-11-29 09:34:17 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDigitizer ------
//
//           by F.Longo, R.Giannitrapani & G.Santin (24 oct 2001)
//
// ************************************************************

#ifndef GammaRayTelDigitizer_h
#define GammaRayTelDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "GammaRayTelDigi.hh"
#include "globals.hh"
//#include "g4std/vector"

class GammaRayTelDigitizerMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelDigitizer : public G4VDigitizerModule
{
public:
  
  GammaRayTelDigitizer(G4String name);
  ~GammaRayTelDigitizer();
  
  void Digitize();
  void SetThreshold(G4double val) { Energy_Threshold = val;}
  
private:
  
  GammaRayTelDigitsCollection*  DigitsCollection;
  G4double Energy_Threshold;
  GammaRayTelDigitizerMessenger* digiMessenger;

};

#endif








