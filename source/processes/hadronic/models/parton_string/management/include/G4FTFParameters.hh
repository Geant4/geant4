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
#ifndef G4FTFCrossSection_h
#define G4FTFCrossSection_h 1
//
// $Id: G4FTFParameters.hh,v 1.1 2008-03-31 14:50:12 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Proton.hh"
class G4FTFparameters
{

  public:
  	G4FTFparameters(const G4ParticleDefinition * , G4double );

	~G4FTFparameters();

        void SethNcmsEnergy(const G4double s);
        void SetTotalCrossSection(const G4double Xtotal);
        void SetElastisCrossSection(const G4double Xelastic);
        void SetInelasticCrossSection(const G4double Xinelastic);

        void SetRadiusOfHNinteractions2(const G4double Radius2);

        void SetSlope(const G4double Slope);
        void SetAvaragePt2ofElasticScattering(const G4double aPt2);
        void SetGamma0(const G4double Gamma0);

	G4double GetTotalCrossSection();
	G4double GetElasticCrossSection();
	G4double GetInelasticCrossSection();
        G4double GetProbabilityOfElasticScatt();

        G4double GetSlope();
        G4double GetAvaragePt2ofElasticScattering();

	G4double GetProbabilityOfInteraction(const G4double impactsquare);
	G4double GetInelasticProbability(const G4double impactsquare);

//  private: 

	G4FTFparameters();

        G4double GammaElastic(const G4double impactsquare) {return (FTFGamma0 * std::exp(-FTFSlope * impactsquare));};
	
        G4double FTFhNcmsEnergy;                // Uzhi Initial hN CMS energy
        G4double FTFXtotal;                     // Total X in mb
        G4double FTFXelastic;                   // Elastic X in mb
        G4double FTFXinelastic;                 // Inelastic X in mb
        G4double ProbabilityOfElasticScatt;     // Xel/Xtot
        G4double RadiusOfHNinteractions2;       // Xtot/pi, in fn^2
        G4double FTFSlope;                      // in fm^-1
        G4double AvaragePt2ofElasticScattering; // in MeV^2
        G4double FTFGamma0;

};

inline  void G4FTFparameters::SethNcmsEnergy(const G4double s)
             {FTFhNcmsEnergy = s;}

inline  void G4FTFparameters::SetTotalCrossSection(const G4double Xtotal)
             {FTFXtotal = Xtotal;}

inline  void G4FTFparameters::SetElastisCrossSection(const G4double Xelastic)
             {FTFXelastic = Xelastic;}

inline  void G4FTFparameters::SetInelasticCrossSection(const G4double Xinelastic)
             {FTFXinelastic = Xinelastic;}

inline  void G4FTFparameters::SetRadiusOfHNinteractions2(const G4double Radius2)
             {RadiusOfHNinteractions2 = Radius2;}

inline  void G4FTFparameters::SetSlope(const G4double Slope)
             {FTFSlope = 12.84/Slope;} // Slope is in GeV^-2, FTFSlope in fm^-2

inline  void G4FTFparameters::SetAvaragePt2ofElasticScattering(const G4double aPt2)
                 {AvaragePt2ofElasticScattering = aPt2;}

inline  void G4FTFparameters::SetGamma0(const G4double Gamma0)
             {FTFGamma0 = Gamma0;}

inline  G4double G4FTFparameters::GetTotalCrossSection()     {return FTFXtotal;}
inline  G4double G4FTFparameters::GetElasticCrossSection()   {return FTFXelastic;}
inline  G4double G4FTFparameters::GetInelasticCrossSection() {return FTFXinelastic;}

inline  G4double G4FTFparameters::GetProbabilityOfElasticScatt()
                 { 
                  if(FTFXtotal==0.) {return 0.;}
                  else              {return FTFXelastic/FTFXtotal;}
                 } 

inline  G4double G4FTFparameters::GetSlope()                 {return FTFSlope;}

inline  G4double G4FTFparameters::GetAvaragePt2ofElasticScattering()
                 {return AvaragePt2ofElasticScattering;}

inline  G4double G4FTFparameters::GetProbabilityOfInteraction(const G4double impactsquare)
                 {
                  if(RadiusOfHNinteractions2 > impactsquare) {return 1.;}
                  else                                       {return 0.;}
                 } 

inline  G4double G4FTFparameters::GetInelasticProbability( const G4double impactsquare)
        {
         G4double Gamma = GammaElastic(impactsquare);
         return 2 * Gamma - Gamma *Gamma;
        }

#endif
