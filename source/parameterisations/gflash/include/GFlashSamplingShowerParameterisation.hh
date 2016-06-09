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
// $Id: GFlashSamplingShowerParameterisation.hh,v 1.3 2005/11/30 19:29:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashSamplingShowerParameterisation
//
//  Class description:
//
//  GFlash concrete sampling shower parameterisation

// Author: Joanna Weng - 02.2004
//---------------------------------------------------------------
#ifndef GFlashSamplingShowerParameterisation_h
#define GFlashSamplingShowerParameterisation_h 1

#include "globals.hh"
#include "GFlashSamplingShowerTuning.hh"
#include "GVFlashShowerParameterisation.hh"

class G4Material;

class GFlashSamplingShowerParameterisation
  : public GVFlashShowerParameterisation
{
  public:

    GFlashSamplingShowerParameterisation(G4Material* aMat1, G4Material* aMat2,
                                         G4double d1, G4double d2,
    GFlashSamplingShowerTuning * aPar = 0);
    ~GFlashSamplingShowerParameterisation();

    void ComputeRadialParameters(G4double y, G4double Tau);
    void GenerateLongitudinalProfile(G4double Energy); 
    void ComputeZAX0EFFetc();

    G4double IntegrateEneLongitudinal(G4double LongitudinalStep);
    G4double IntegrateNspLongitudinal(G4double LongitudinalStep);
    G4double ComputeTau(G4double LongitudinalPosition);
    void SetMaterial(G4Material *mat1, G4Material *mat2);
    G4double GeneratePhi();
    G4double GenerateRadius(G4int ispot, G4double Energy,
    G4double LongitudinalPosition);
    G4double GenerateExponential(G4double Energy);

    inline G4double GetAveR99() {return (3.5 * Rmeff);}
    inline G4double GetAveR90() {return (1.5 * Rmeff);} //ok
    //
    inline G4double GetAveTmx() {return (X0eff*std::exp(AveLogTmax));}
    inline G4double GetAveT99() {return (X0eff*AveLogTmax/(AveLogAlpha-1.00));}
    inline G4double GetAveT90() {return (2.5* X0eff* std::exp( AveLogTmax));}
    // 
    inline G4double GetNspot()  {return NSpot;}
    inline G4double GetX0()     {return X0eff;}  
    inline G4double GetEc()     {return Eceff;} 
    inline G4double GetRm()     {return Rmeff;} 

    G4double ApplySampling(const G4double DEne, const G4double Energy);

  private:

    // medium related quantities
    //
    G4Material *material1, *material2 ;
    G4double  density1, A1, Z1, X01, Ec1, Rm1, d1;
    G4double  density2, A2, Z2, X02, Ec2, Rm2, d2;
    G4double  Aeff, Rhoeff, X0eff, Eceff, Rmeff, Fs, ehat, Zeff;

    // Resolution
    //
    G4double ConstantResolution; 
    G4double NoiseResolution;   
    G4double SamplingResolution;
     
    // parametrization parameters
    //
    GFlashSamplingShowerTuning * thePar;

    // Cashed parameters:  
    // Longitudinal Coefficients for a homogenious calo
    //
    G4double ParAveT1, ParAveT2;
    G4double ParAveA1,ParAveA2, ParAveA3;
    G4double ParSigLogT1,ParSigLogT2;
    G4double ParSigLogA1,ParSigLogA2;
    G4double ParRho1,ParRho2;

    //Cashed parameters:  
    // Longitudinal Coefficients for a sampling calo
    //
    G4double ParsAveT1, ParsAveT2;
    G4double ParsAveA1,ParsAveA2;
    G4double ParsSigLogT1,ParsSigLogT2;
    G4double ParsSigLogA1,ParsSigLogA2;
    G4double ParsRho1,ParsRho2;
    void ComputeLongitudinalParameters(G4double y);
    void GenerateEnergyProfile(G4double y);
    void GenerateNSpotProfile(G4double y);

    // Radial Coefficients homo
    //
    G4double ParRC1,ParRC2,ParRC3,ParRC4;
    G4double ParWC1,ParWC2,ParWC3;
    G4double ParWC4,ParWC5,ParWC6;
    G4double ParRT1,ParRT2,ParRT3,ParRT4;
    G4double ParRT5,ParRT6;

    // Radial Coefficients sampling
    //
    G4double ParsRC1,ParsRC2;
    G4double ParsWC1,ParsWC2;   
    G4double ParsRT1,ParsRT2;

    // Spot multiplicity Coefficients
    //
    G4double ParsSpotT1,ParsSpotT2,ParsSpotA1, ParsSpotA2;
    G4double ParsSpotN1,ParsSpotN2;

    // PARAMETRISATION variables (Energy & position dependent)
    // Longitudinal 
    // homogeneous
    //
    G4double AveLogAlphah,AveLogTmaxh;
    G4double SigmaLogAlphah,SigmaLogTmaxh;
    G4double Rhoh;
    G4double Alphah,Tmaxh,Betah;  

    // PARAMETRISATION variables (Energy & position dependent)
    // Longitudinal 
    // sampling
    //
    G4double AveLogAlpha,AveLogTmax;
    G4double SigmaLogAlpha,SigmaLogTmax;
    G4double Rho;
    G4double Alpha,Tmax,Beta;  

    // Multiplicity
    //
    G4double NSpot,AlphaNSpot,TNSpot,BetaNSpot;

    //Radial
    //
    G4double RadiusCore, WeightCore,RadiusTail; 
};

#endif
