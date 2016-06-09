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
// Created by J. Weng  02.2004
#ifndef GFlashHomoShowerParamterisation_h
#define GFlashHomoShowerParamterisation_h 1

#include "globals.hh"
#include "GVFlashHomoShowerTuning.hh"

class G4Material;
class  GFlashHomoShowerParamterisation
{
	public:
	GFlashHomoShowerParamterisation(G4Material * aMat, GVFlashHomoShowerTuning * aPar = 0);
	~GFlashHomoShowerParamterisation();
	
	void ComputeRadialParameters(G4double y, G4double Tau);
	void GenerateLongitudinalProfile(G4double Energy); 
	void ComputeZAX0EFFetc();
	
	G4double IntegrateEneLongitudinal(G4double LongitudinalStep);
	G4double IntegrateNspLongitudinal(G4double LongitudinalStep);
	G4double ComputeTau(G4double LongitudinalPosition);
	
	G4double GeneratePhi();
	G4double GenerateRadius(G4int ispot, G4double Energy, G4double LongitudinalPosition);
	G4double GenerateExponential(G4double Energy);
	
	//Gets
	//R
	inline   G4double  GetAveR99() {return (3.5 * Rm);}
	inline   G4double  GetAveR90() {return (1.5 * Rm);} //ok
	//T
	inline   G4double  GetAveTmx() {return (X0 *  std::exp(AveLogTmaxh));}  //exp ?
	inline   G4double  GetAveT99() {return (X0 *  AveLogTmaxh/(AveLogAlphah-1.00));}
	inline   G4double  GetAveT90() {return (2.5 * X0 * std::exp( AveLogTmaxh) );}
	// 
	inline   G4double GetNspot(){ return NSpot;}
	G4double GetEffZ(const G4Material * material);
	G4double GetEffA(const G4Material * material);   
	G4double gam(G4double x, G4double a) const; // @@@@ gamma function
	
	private:
	// medium related quantities
	G4double  density, A, Z, X0, Ec, Rm;
	G4Material *material;
	
	public:
	inline   G4double GetX0(){return X0;}  
	inline   G4double GetEc(){return Ec;} 
	inline   G4double GetRm(){return Rm;} 
	void SetMaterial (G4Material *mat);
	void PrintMaterial();
	
	private:
	//Resolution
	G4double ConstantResolution; 
	G4double NoiseResolution;   
	G4double SamplingResolution;
	
	// parametrization parameters
	GVFlashHomoShowerTuning * thePar;
	
	// Cashed parameters:  
	// Longitudinal Coefficients for a homogenious calo
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
	
	//Spot multiplicity Coefficients
	G4double ParSpotT1,ParSpotT2,ParSpotA1, ParSpotA2;
	G4double ParSpotN1,ParSpotN2;
	
	//PARAMETRISATION variables (Energy & position dependent)
	//Longitudinal 
	//homogenous
	G4double AveLogAlphah,AveLogTmaxh;
	G4double SigmaLogAlphah,SigmaLogTmaxh;
	G4double Rhoh;
	G4double Alphah,Tmaxh,Betah;  
	
	//MMultiplicity
	G4double NSpot,AlphaNSpot,TNSpot,BetaNSpot;
	
	//Radial
	G4double RadiusCore, WeightCore,RadiusTail; 
};

#endif

