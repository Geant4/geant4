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
// $Id: G4PAIySection.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// 
// G4PAIySection.hh -- header file
//
//
// Preparation of ionizing collision cross section according to Photo Absorption 
// Ionization (PAI) model for simulation of ionization energy losses in very thin
// absorbers. Author: Vladimir.Grichine@cern.ch
//
// History:
//
// 01.10.07, V.Ivanchenko create using V.Grichine G4PAIxSection class
// 21.11.10, V.Grichine   fVerbose and SetVerbose added
// 28.10.11, V.Ivanchenko Migration of exceptions to the new design 

#ifndef G4PAIYSECTION_HH
#define G4PAIYSECTION_HH

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4SandiaTable.hh"

class G4PAIySection
{
public:

  explicit G4PAIySection();
	  
  ~G4PAIySection();

  void Initialize(const G4Material* material, G4double maxEnergyTransfer, 
		  G4double betaGammaSq, G4SandiaTable*);

  void ComputeLowEnergyCof(const G4Material* material);

  void InitPAI();

  void NormShift( G4double betaGammaSq );
  
  void SplainPAI( G4double betaGammaSq );
	  	  
  // Physical methods
  G4double RutherfordIntegral( G4int intervalNumber,
			       G4double limitLow,
			       G4double limitHigh     );

  G4double ImPartDielectricConst( G4int intervalNumber,
				  G4double energy        );

  G4double RePartDielectricConst(G4double energy);

  G4double DifPAIySection( G4int intervalNumber,
			   G4double betaGammaSq    );

  G4double PAIdNdxCerenkov( G4int intervalNumber,
			    G4double betaGammaSq    );

  G4double PAIdNdxPlasmon( G4int intervalNumber,
			   G4double betaGammaSq    );

  void     IntegralPAIySection();
  void     IntegralCerenkov();
  void     IntegralPlasmon();

  G4double SumOverInterval(G4int intervalNumber);
  G4double SumOverIntervaldEdx(G4int intervalNumber);
  G4double SumOverInterCerenkov(G4int intervalNumber);
  G4double SumOverInterPlasmon(G4int intervalNumber);

  G4double SumOverBorder( G4int intervalNumber,
			  G4double energy          );
  G4double SumOverBorderdEdx( G4int intervalNumber,
			      G4double energy          );
  G4double SumOverBordCerenkov( G4int intervalNumber,
				G4double energy          );
  G4double SumOverBordPlasmon( G4int intervalNumber,
			       G4double energy          );

  G4double GetStepEnergyLoss( G4double step );
  G4double GetStepCerenkovLoss( G4double step );
  G4double GetStepPlasmonLoss( G4double step );
	 
  // Inline access functions

  inline G4int GetNumberOfGammas() const { return fNumberOfGammas; }
	  
  inline G4int GetSplineSize() const { return fSplineNumber; }
	  
  inline G4int GetIntervalNumber() const { return fIntervalNumber; }

  inline G4double GetEnergyInterval(G4int i){ return fEnergyInterval[i]; } 

  inline G4double GetDifPAIySection(G4int i){ return fDifPAIySection[i]; } 
  inline G4double GetPAIdNdxCrenkov(G4int i){ return fdNdxCerenkov[i]; } 
  inline G4double GetPAIdNdxPlasmon(G4int i){ return fdNdxPlasmon[i]; } 
	  
  inline G4double GetMeanEnergyLoss() const {return fIntegralPAIySection[0]; }
  inline G4double GetMeanCerenkovLoss() const {return fIntegralCerenkov[0]; }
  inline G4double GetMeanPlasmonLoss() const {return fIntegralPlasmon[0]; }

  inline G4double GetNormalizationCof() const { return fNormalizationCof; }
          
  inline G4double GetPAItable(G4int i,G4int j) const;

  inline G4double GetLorentzFactor(G4int i) const;
	  	  
  inline G4double GetSplineEnergy(G4int i) const;
	  
  inline G4double GetIntegralPAIySection(G4int i) const;
  inline G4double GetIntegralPAIdEdx(G4int i) const;
  inline G4double GetIntegralCerenkov(G4int i) const;
  inline G4double GetIntegralPlasmon(G4int i) const;

  inline void SetVerbose(G4int v) { fVerbose = v; };

private :

  void CallError(G4int i, const G4String& methodName) const;

  // Local class constants
 
  static const G4double fDelta; // energy shift from interval border = 0.001
  static const G4double fError; // error in lin-log approximation = 0.005

  static G4int fNumberOfGammas; // = 111;
  static const G4double fLorentzFactor[112];  //  static gamma array

  static
  const G4int fRefGammaNumber; // The number of gamma for creation of spline (15)

  G4int    fIntervalNumber ;   //  The number of energy intervals
  G4double fNormalizationCof;  // Normalization cof for PhotoAbsorptionXsection

  G4double betaBohr;
  G4double betaBohr4;

  G4double fDensity;            // Current density
  G4double fElectronDensity;    // Current electron (number) density
  G4double fLowEnergyCof;       // Correction cof for low energy region
  G4int    fSplineNumber;       // Current size of spline
  G4int    fVerbose;            // verbose flag

  G4SandiaTable*  fSandia;

  G4DataVector fEnergyInterval;
  G4DataVector fA1; 
  G4DataVector fA2;
  G4DataVector fA3; 
  G4DataVector fA4;

  static
  const G4int  fMaxSplineSize; // Max size of output splain arrays = 500

  G4DataVector fSplineEnergy;          // energy points of splain
  G4DataVector fRePartDielectricConst; // Real part of dielectric const
  G4DataVector fImPartDielectricConst; // Imaginary part of dielectric const
  G4DataVector fIntegralTerm;          // Integral term in PAI cross section
  G4DataVector fDifPAIySection;        // Differential PAI cross section
  G4DataVector fdNdxCerenkov;          // dNdx of Cerenkov collisions
  G4DataVector fdNdxPlasmon;           // dNdx of Plasmon collisions

  G4DataVector fIntegralPAIySection;   // Integral PAI cross section  ?
  G4DataVector fIntegralPAIdEdx;       // Integral PAI dEdx  ?
  G4DataVector fIntegralCerenkov;      // Integral Cerenkov N>omega  ?
  G4DataVector fIntegralPlasmon;       // Integral Plasmon N>omega  ?

  G4double     fPAItable[500][112];    // Output array
};

inline G4double G4PAIySection::GetPAItable(G4int i, G4int j) const
{
   return fPAItable[i][j];
}

inline G4double G4PAIySection::GetLorentzFactor(G4int j) const
{
   return fLorentzFactor[j];
}

inline G4double G4PAIySection::GetSplineEnergy(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetSplineEnergy"); }
  return fSplineEnergy[i];
}
	  
inline G4double G4PAIySection::GetIntegralPAIySection(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPAIySection"); }
  return fIntegralPAIySection[i];
}

inline G4double G4PAIySection::GetIntegralPAIdEdx(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPAIdEdx"); }
  return fIntegralPAIdEdx[i];
}

inline G4double G4PAIySection::GetIntegralCerenkov(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralCerenkov"); }
  return fIntegralCerenkov[i];
}

inline G4double G4PAIySection::GetIntegralPlasmon(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPlasmon"); }
  return fIntegralPlasmon[i];
}

#endif   

// -----------------   end of G4PAIySection header file    -------------------
