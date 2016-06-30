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
// $Id: G4PAIxSection.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// 
// G4PAIxSection.hh -- header file
//
// GEANT 4 class header file --- Copyright CERN 1995
// CERB Geneva Switzerland
//
// for information related to this code, please, contact
// CERN, CN Division, ASD Group
//
// Preparation of ionizing collision cross section according to Photo Absorption 
// Ionization (PAI) model for simulation of ionization energy losses in very thin
// absorbers. Author: Vladimir.Grichine@cern.ch
//
// History:
//
// 28.10.11, V. Ivanchenko: Migration of exceptions to the new design 
// 19.10.03, V. Grichine: Integral dEdx was added for G4PAIModel class 
// 13.05.03, V. Grichine: Numerical instability was fixed in SumOverInterval/Border 
//                        functions
// 10.02.02, V. Grichine: New functions and arrays/gets for Cerenkov and 
//                        plasmon collisions dN/dx
// 27.10.99, V. Grichine: Bug fixed in constructors, 3rd constructor and 
//                        GetStepEnergyLoss(step) were added, fDelta = 0.005
// 30.11.97, V. Grichine: 2nd version 
// 11.06.97, V. Grichine: 1st version

#ifndef G4PAIXSECTION_HH
#define G4PAIXSECTION_HH

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

#include"G4SandiaTable.hh"

class G4MaterialCutsCouple;
class G4Sandiatable;


class G4PAIxSection
{
public:
	  // Constructors
  G4PAIxSection();
  G4PAIxSection( G4MaterialCutsCouple* matCC);
	  
  G4PAIxSection( G4int materialIndex,
	                 G4double maxEnergyTransfer   );
	  
  G4PAIxSection( G4int materialIndex,           // for proton loss table
		         G4double maxEnergyTransfer,
		         G4double betaGammaSq ,
                         G4double** photoAbsCof, G4int intNumber         );

  G4PAIxSection( G4int materialIndex,           // test constructor
		         G4double maxEnergyTransfer,
		         G4double betaGammaSq          );
	  
  ~G4PAIxSection();

  void Initialize(const G4Material* material, G4double maxEnergyTransfer, 
		  G4double betaGammaSq, G4SandiaTable*);
	  
  // General control functions

  void     ComputeLowEnergyCof(const G4Material* material);
  void     ComputeLowEnergyCof();
          
  void InitPAI();

  void NormShift( G4double betaGammaSq );

  void SplainPAI( G4double betaGammaSq );
	  	  
  // Physical methods
          
  G4double RutherfordIntegral( G4int intervalNumber,
			       G4double limitLow,
			       G4double limitHigh     );

  G4double ImPartDielectricConst( G4int intervalNumber,
				  G4double energy        );

  G4double GetPhotonRange( G4double energy );
  G4double GetElectronRange( G4double energy );

  G4double RePartDielectricConst(G4double energy);

  G4double DifPAIxSection( G4int intervalNumber,
			   G4double betaGammaSq    );

  G4double PAIdNdxCerenkov( G4int intervalNumber,
			    G4double betaGammaSq    );
  G4double PAIdNdxMM( G4int intervalNumber,
		      G4double betaGammaSq    );

  G4double PAIdNdxPlasmon( G4int intervalNumber,
			   G4double betaGammaSq    );

  G4double PAIdNdxResonance( G4int intervalNumber,
			     G4double betaGammaSq    );


  void     IntegralPAIxSection();
  void     IntegralCerenkov();
  void     IntegralMM();
  void     IntegralPlasmon();
  void     IntegralResonance();

  G4double SumOverInterval(G4int intervalNumber);
  G4double SumOverIntervaldEdx(G4int intervalNumber);
  G4double SumOverInterCerenkov(G4int intervalNumber);
  G4double SumOverInterMM(G4int intervalNumber);
  G4double SumOverInterPlasmon(G4int intervalNumber);
  G4double SumOverInterResonance(G4int intervalNumber);

  G4double SumOverBorder( G4int intervalNumber,
			  G4double energy          );
  G4double SumOverBorderdEdx( G4int intervalNumber,
			      G4double energy          );
  G4double SumOverBordCerenkov( G4int intervalNumber,
				G4double energy          );
  G4double SumOverBordMM( G4int intervalNumber,
			  G4double energy          );
  G4double SumOverBordPlasmon( G4int intervalNumber,
			       G4double energy          );
  G4double SumOverBordResonance( G4int intervalNumber,
				 G4double energy          );

  G4double GetStepEnergyLoss( G4double step );
  G4double GetStepCerenkovLoss( G4double step );
  G4double GetStepMMLoss( G4double step );
  G4double GetStepPlasmonLoss( G4double step );
  G4double GetStepResonanceLoss( G4double step );
	 
  G4double GetEnergyTransfer();
  G4double GetCerenkovEnergyTransfer();
  G4double GetMMEnergyTransfer();
  G4double GetPlasmonEnergyTransfer();
  G4double GetResonanceEnergyTransfer();
  G4double GetRutherfordEnergyTransfer();
	 
  // Inline access functions

  G4int GetNumberOfGammas() const { return fNumberOfGammas; }
	  
  G4int GetSplineSize() const { return fSplineNumber; }
	  
  G4int GetIntervalNumber() const { return fIntervalNumber; }

  G4double GetEnergyInterval(G4int i){ return fEnergyInterval[i]; } 

  G4double GetDifPAIxSection(G4int i){ return fDifPAIxSection[i]; } 
  G4double GetPAIdNdxCerenkov(G4int i){ return fdNdxCerenkov[i]; } 
  G4double GetPAIdNdxMM(G4int i){ return fdNdxMM[i]; } 
  G4double GetPAIdNdxPlasmon(G4int i){ return fdNdxPlasmon[i]; } 
  G4double GetPAIdNdxResonance(G4int i){ return fdNdxResonance[i]; } 
	  
  G4double GetMeanEnergyLoss() const {return fIntegralPAIxSection[0]; }
  G4double GetMeanCerenkovLoss() const {return fIntegralCerenkov[0]; }
  G4double GetMeanMMLoss() const {return fIntegralMM[0]; }
  G4double GetMeanPlasmonLoss() const {return fIntegralPlasmon[0]; }
  G4double GetMeanResonanceLoss() const {return fIntegralResonance[0]; }

  G4double GetNormalizationCof() const { return fNormalizationCof; }

  G4double GetLowEnergyCof() const { return fLowEnergyCof; }

  void SetVerbose(G4int v){fVerbose=v;};
          
  inline G4double GetPAItable(G4int i,G4int j) const;
  
  inline G4double GetLorentzFactor(G4int i) const;
	  	  
  inline G4double GetSplineEnergy(G4int i) const;
	  
  inline G4double GetIntegralPAIxSection(G4int i) const;
  inline G4double GetIntegralPAIdEdx(G4int i) const;
  inline G4double GetIntegralCerenkov(G4int i) const;
  inline G4double GetIntegralMM(G4int i) const;
  inline G4double GetIntegralPlasmon(G4int i) const;
  inline G4double GetIntegralResonance(G4int i) const;

private :

  void CallError(G4int i, const G4String& methodName) const;

  G4PAIxSection & operator=(const G4PAIxSection &right) = delete;
  G4PAIxSection(const G4PAIxSection&) = delete;

  // Local class constants
 
  static const G4double fDelta; // energy shift from interval border = 0.001
  static const G4double fError; // error in lin-log approximation = 0.005

  static G4int fNumberOfGammas;         // = 111;
  static const G4double fLorentzFactor[112];  //  static gamma array

  static 
  const G4int fRefGammaNumber; //The number of gamma for creation of spline (15)

  G4int    fIntervalNumber ;    //  The number of energy intervals
  G4double fNormalizationCof;   // Normalization cof for PhotoAbsorptionXsection

  // G4double fBetaGammaSq;        // (beta*gamma)^2

  G4int fMaterialIndex;  // current material index
  G4double fDensity;            // Current density
  G4double fElectronDensity;    // Current electron (number) density
  G4double fLowEnergyCof;    // Correction cof for low energy region
  G4int    fSplineNumber;       // Current size of spline
  G4int    fVerbose;       // verbose flag

  // Arrays of Sandia coefficients

  G4OrderedTable* fMatSandiaMatrix;

  G4SandiaTable*  fSandia;

  G4DataVector fEnergyInterval;
  G4DataVector fA1; 
  G4DataVector fA2;
  G4DataVector fA3; 
  G4DataVector fA4;

  static
  const G4int  fMaxSplineSize ; // Max size of output splain arrays = 500

  G4DataVector          fSplineEnergy;   // energy points of splain
  G4DataVector   fRePartDielectricConst;   // Real part of dielectric const
  G4DataVector   fImPartDielectricConst;   // Imaginary part of dielectric const
  G4DataVector            fIntegralTerm;   // Integral term in PAI cross section
  G4DataVector          fDifPAIxSection;   // Differential PAI cross section
  G4DataVector            fdNdxCerenkov;   // dNdx of Cerenkov collisions
  G4DataVector            fdNdxPlasmon;    // dNdx of Plasmon collisions
  G4DataVector            fdNdxMM;    // dNdx of MM-Cerenkov collisions
  G4DataVector            fdNdxResonance;    // dNdx of Resonance collisions

  G4DataVector     fIntegralPAIxSection;   // Integral PAI cross section  ?
  G4DataVector     fIntegralPAIdEdx;       // Integral PAI dEdx  ?
  G4DataVector     fIntegralCerenkov;      // Integral Cerenkov N>omega  ?
  G4DataVector     fIntegralPlasmon;       // Integral Plasmon N>omega  ?
  G4DataVector     fIntegralMM;       // Integral MM N>omega  ?
  G4DataVector     fIntegralResonance;       // Integral resonance N>omega  ?

  G4double fPAItable[500][112]; // Output array

};    

////////////////  Inline methods //////////////////////////////////
//

inline G4double G4PAIxSection::GetPAItable(G4int i, G4int j) const
{
   return fPAItable[i][j];
}

inline G4double G4PAIxSection::GetLorentzFactor(G4int j) const
{
   return fLorentzFactor[j];
}

inline G4double G4PAIxSection::GetSplineEnergy(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetSplineEnergy"); }
  return fSplineEnergy[i];
}
	  
inline G4double G4PAIxSection::GetIntegralPAIxSection(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPAIxSection"); }
  return fIntegralPAIxSection[i];
}

inline G4double G4PAIxSection::GetIntegralPAIdEdx(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPAIdEdx"); }
  return fIntegralPAIdEdx[i];
}

inline G4double G4PAIxSection::GetIntegralCerenkov(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralCerenkov"); }
  return fIntegralCerenkov[i];
}

inline G4double G4PAIxSection::GetIntegralMM(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralMM"); }
  return fIntegralMM[i];
}

inline G4double G4PAIxSection::GetIntegralPlasmon(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralPlasmon"); }
  return fIntegralPlasmon[i];
}

inline G4double G4PAIxSection::GetIntegralResonance(G4int i) const 
{
  if(i < 1 || i > fSplineNumber) { CallError(i, "GetIntegralResonance"); }
  return fIntegralResonance[i];
}

#endif   

// -----------------   end of G4PAIxSection header file    -------------------
