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
// $Id: G4FTFCrossSection.hh,v 1.2 2007/04/24 10:37:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
#include "G4Proton.hh"
class G4FTFCrossSection
{

  public:
  	G4FTFCrossSection(const G4ParticleDefinition * , G4double );

	~G4FTFCrossSection();

        void SethNcmsEnergy(const G4double s);
        void SetTotalCrossSection(const G4double Xtotal);
        void SetElastisCrossSection(const G4double Xelastic);
        void SetInelasticCrossSection(const G4double Xinelastic);
        void SetSlope(const G4double Slope);
        void SetGamma0(const G4double Gamma0);

	G4double GetTotalCrossSection();
	G4double GetElasticCrossSection();
	G4double GetInelasticCrossSection();

        G4double GetSlope();

	G4double GetInelasticProbability(const G4double impactsquare);

//  private: 

	G4FTFCrossSection();

        G4double GammaElastic(const G4double impactsquare) {return (FTFGamma0 * std::exp(-FTFSlope * impactsquare));};
	
        G4double FTFhNcmsEnergy;  // Uzhi Initial hN CMS energy
        G4double FTFXtotal;
        G4double FTFXelastic;
        G4double FTFXinelastic;
        G4double FTFSlope;
        G4double FTFGamma0;

};

inline  void G4FTFCrossSection::SethNcmsEnergy(const G4double s)                    {FTFhNcmsEnergy = s;}
inline  void G4FTFCrossSection::SetTotalCrossSection(const G4double Xtotal)         {FTFXtotal = Xtotal;}
inline  void G4FTFCrossSection::SetElastisCrossSection(const G4double Xelastic)     {FTFXelastic = Xelastic;}
inline  void G4FTFCrossSection::SetInelasticCrossSection(const G4double Xinelastic) {FTFXinelastic = Xinelastic;}
inline  void G4FTFCrossSection::SetSlope(const G4double Slope)                      {FTFSlope = 12.84/Slope;}
inline  void G4FTFCrossSection::SetGamma0(const G4double Gamma0)                    {FTFGamma0 = Gamma0;}

inline  G4double G4FTFCrossSection::GetTotalCrossSection()                          {return FTFXtotal;}
inline  G4double G4FTFCrossSection::GetElasticCrossSection()                        {return FTFXelastic;}
inline  G4double G4FTFCrossSection::GetInelasticCrossSection()                      {return FTFXinelastic;}

inline  G4double G4FTFCrossSection::GetSlope()                                      {return FTFSlope;}

inline  G4double G4FTFCrossSection::GetInelasticProbability( const G4double impactsquare)
        {
         G4double Gamma = GammaElastic(impactsquare);
         return 2 * Gamma - Gamma *Gamma;
        }

#endif
