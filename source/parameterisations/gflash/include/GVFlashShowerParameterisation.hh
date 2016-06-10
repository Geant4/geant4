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
// $Id: GVFlashShowerParameterisation.hh 93079 2015-10-02 14:43:41Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GVFlashShowerParameterisation
//
//  Class description:
//
//  Base class for GFlash shower parameterisation.

// Author: Joanna Weng - 11.2005
//---------------------------------------------------------------
#ifndef GVFlashShowerParameterisation_h
#define GVFlashShowerParameterisation_h 1

#include "globals.hh"
#include "GVFlashHomoShowerTuning.hh"

class MyGamma;
class G4Material;


class GVFlashShowerParameterisation
{
public:  // with description
  
  GVFlashShowerParameterisation();
  virtual ~GVFlashShowerParameterisation();
  
  virtual void ComputeRadialParameters(G4double y, G4double Tau)       = 0;
  virtual void GenerateLongitudinalProfile(G4double Energy)            = 0; 
  virtual G4double IntegrateEneLongitudinal(G4double LongitudinalStep) = 0;
  virtual G4double IntegrateNspLongitudinal(G4double LongitudinalStep) = 0;
  virtual G4double ComputeTau(G4double LongitudinalPosition)           = 0;
  virtual G4double GenerateRadius(G4int ispot, G4double Energy,
				  G4double LongitudinalPosition)       = 0; 
  virtual void ComputeLongitudinalParameters(G4double y)               = 0;
  virtual void GenerateEnergyProfile(G4double y)                       = 0;
  virtual void GenerateNSpotProfile(G4double y)                        = 0;
  virtual G4double GenerateExponential(G4double Energy)                = 0;
  
  virtual G4double GetAveR99() = 0;
  virtual G4double GetAveR90() = 0;
  
  virtual G4double GetAveTmx() = 0;
  virtual G4double GetAveT99() = 0; 
  virtual G4double GetAveT90() = 0; 
  
  virtual G4double GetNspot()  = 0;
  virtual G4double GetX0()     = 0;
  virtual G4double GetEc()     = 0;
  virtual G4double GetRm()     = 0;
  
  G4double GeneratePhi();
  G4double GetEffZ(const G4Material * material);
  G4double GetEffA(const G4Material * material);   
  G4double gam(G4double x, G4double a) const; // @@@@ gamma function
  void PrintMaterial(const G4Material * mat);
  
protected:
  
  GVFlashHomoShowerTuning * thePar;
  // Parameterisation parameters
  G4double  density, A, Z, X0, Ec, Rm;
  // Medium related quantities
  G4double NSpot;
  
private:
  MyGamma* fGamma;
  
};

#endif
