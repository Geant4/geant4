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
// $Id: GVFlashShowerParameterisation.hh,v 1.1 2005/11/30 19:29:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
};

#endif
