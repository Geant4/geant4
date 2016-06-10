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
// $Id: GFlashHomoShowerParameterisation.hh 68057 2013-03-13 14:46:00Z gcosmo $
//
//
// ---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashHomoShowerParameterisation
//
//  Class description:
//
//  GFlash homogeneous shower parameterisation.

// Author: Joanna Weng - 02.2004
// ---------------------------------------------------------------
#ifndef GFlashHomoShowerParameterisation_h
#define GFlashHomoShowerParameterisation_h 1

#include "globals.hh"
#include "GVFlashHomoShowerTuning.hh"
#include "GVFlashShowerParameterisation.hh"

class G4Material;

class GFlashHomoShowerParameterisation : public GVFlashShowerParameterisation
{
  public:  // with description
  
    GFlashHomoShowerParameterisation(G4Material * aMat,
                                     GVFlashHomoShowerTuning * aPar = 0);
    ~GFlashHomoShowerParameterisation();

    void ComputeRadialParameters(G4double y, G4double Tau);
    void GenerateLongitudinalProfile(G4double Energy); 
    void ComputeZAX0EFFetc();

    G4double IntegrateEneLongitudinal(G4double LongitudinalStep);
    G4double IntegrateNspLongitudinal(G4double LongitudinalStep);
    G4double ComputeTau(G4double LongitudinalPosition);

    G4double GeneratePhi();
    G4double GenerateRadius(G4int ispot, G4double Energy,
    G4double LongitudinalPosition);
    G4double GenerateExponential(G4double Energy);
    void SetMaterial(G4Material *mat);

    inline G4double GetAveR99() {return (3.5 * Rm);}
    inline G4double GetAveR90() {return (1.5 * Rm);} //ok

    inline G4double GetAveTmx() {return (X0 * std::exp(AveLogTmaxh));}
    inline G4double GetAveT99() {return (X0 * AveLogTmaxh/(AveLogAlphah-1.00));}
    inline G4double GetAveT90() {return (2.5* X0*std::exp( AveLogTmaxh) );}

    inline   G4double GetNspot(){ return NSpot;}
    inline   G4double GetX0(){return X0;}  
    inline   G4double GetEc(){return Ec;} 
    inline   G4double GetRm(){return Rm;} 

  private:

    G4Material *material;

    //Resolution
    G4double ConstantResolution; 
    G4double NoiseResolution;   
    G4double SamplingResolution;

    // parametrization parameters
    GVFlashHomoShowerTuning * thePar;
    G4bool owning;

    // Cashed parameters:  
    // Longitudinal Coefficients for a homogeneous calo
    G4double ParAveT1;
    G4double ParAveA1,ParAveA2,ParAveA3;
    G4double ParSigLogT1,ParSigLogT2;
    G4double ParSigLogA1,ParSigLogA2;
    G4double ParRho1,ParRho2;

    void ComputeLongitudinalParameters(G4double y);
    void GenerateEnergyProfile(G4double y);
    void GenerateNSpotProfile(G4double y);

    // Radial Coefficients
    G4double ParRC1,ParRC2,ParRC3,ParRC4;
    G4double ParWC1,ParWC2,ParWC3;
    G4double ParWC4,ParWC5,ParWC6;
    G4double ParRT1,ParRT2,ParRT3,ParRT4;
    G4double ParRT5,ParRT6;

    // Spot multiplicity Coefficients
    G4double ParSpotT1,ParSpotT2,ParSpotA1, ParSpotA2;
    G4double ParSpotN1,ParSpotN2;

    // PARAMETRISATION variables (Energy & position dependent)
    // Longitudinal 
    // homogeneous
    G4double AveLogAlphah,AveLogTmaxh;
    G4double SigmaLogAlphah,SigmaLogTmaxh;
    G4double Rhoh;
    G4double Alphah,Tmaxh,Betah;  

    // Multiplicity
    G4double NSpot,AlphaNSpot,TNSpot,BetaNSpot;

    //Radial
    G4double RadiusCore, WeightCore,RadiusTail;
};

#endif

