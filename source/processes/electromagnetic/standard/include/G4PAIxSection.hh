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
// $Id: G4PAIxSection.hh,v 1.9 2002-10-14 17:30:10 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// absorbers
//
// History:
// 1st version 11.06.97, V. Grichine 
// 2nd version 30.11.97, V. Grichine
// 27.10.99, V.Grichine: Bug fixed in constructors, 3rd constructor and 
//                       GetStepEnergyLoss(step) were added, fDelta = 0.005
// 10.02.02, V.Grichine: New functions and arrays/gets for Cerenkov and 
//                       plasmon collisions dN/dx

#ifndef G4PAIXSECTION_HH
#define G4PAIXSECTION_HH

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

#include"G4SandiaTable.hh"

class G4PAIxSection
{
public:
	  // Constructors
	  
	  G4PAIxSection( G4int materialIndex,
	                 G4double maxEnergyTransfer   ) ;
	  
          G4PAIxSection( G4int materialIndex,           // for proton loss table
		         G4double maxEnergyTransfer,
		         G4double betaGammaSq ,
                         G4double** photoAbsCof, G4int intNumber         ) ;

          G4PAIxSection( G4int materialIndex,           // test constructor
		         G4double maxEnergyTransfer,
		         G4double betaGammaSq          ) ;
	  
	  // G4PAIxSection(const G4PAIxSection& right) ;
	  
	  // Destructor
	  
          ~G4PAIxSection() ;
	  
	  // Operators
          // G4PAIxSection& operator=(const G4PAIxSection& right) ;
          // G4int operator==(const G4PAIxSection& right)const ;
          // G4int operator!=(const G4PAIxSection& right)const ;
	  
	  // Methods
	  
	  // General control functions
          
	  void InitPAI() ;

          void NormShift( G4double betaGammaSq ) ;

          void SplainPAI( G4double betaGammaSq ) ;
	  	  
	  // Physical methods
          

          G4double RutherfordIntegral( G4int intervalNumber,
	                               G4double limitLow,
				       G4double limitHigh     ) ;

          G4double ImPartDielectricConst( G4int intervalNumber,
	                                  G4double energy        ) ;

          G4double RePartDielectricConst(G4double energy) ;

          G4double DifPAIxSection( G4int intervalNumber,
	                           G4double betaGammaSq    ) ;

          G4double PAIdNdxCerenkov( G4int intervalNumber,
	                           G4double betaGammaSq    ) ;

          G4double PAIdNdxPlasmon( G4int intervalNumber,
	                           G4double betaGammaSq    ) ;

	  void     IntegralPAIxSection() ;
	  void     IntegralCerenkov() ;
	  void     IntegralPlasmon() ;

          G4double SumOverInterval(G4int intervalNumber) ;
          G4double SumOverInterCerenkov(G4int intervalNumber) ;
          G4double SumOverInterPlasmon(G4int intervalNumber) ;

          G4double SumOverBorder( G4int intervalNumber,
	                          G4double energy          ) ;
          G4double SumOverBordCerenkov( G4int intervalNumber,
	                                G4double energy          ) ;
          G4double SumOverBordPlasmon( G4int intervalNumber,
	                               G4double energy          ) ;

          G4double GetStepEnergyLoss( G4double step ) ;
          G4double GetStepCerenkovLoss( G4double step ) ;
          G4double GetStepPlasmonLoss( G4double step ) ;
	 
	  // Inline access functions

	  G4int GetNumberOfGammas() const { return fNumberOfGammas ; }
	  
          G4int GetSplineSize() const { return fSplineNumber ; }
	  
          G4int GetIntervalNumber() const { return fIntervalNumber ; }

          G4double GetEnergyInterval(G4int i){ return fEnergyInterval[i] ; } 

          G4double GetDifPAIxSection(G4int i){ return fDifPAIxSection[i] ; } 
          G4double GetPAIdNdxCrenkov(G4int i){ return fdNdxCerenkov[i] ; } 
          G4double GetPAIdNdxPlasmon(G4int i){ return fdNdxPlasmon[i] ; } 
	  
	  G4double GetMeanEnergyLoss() const {return fIntegralPAIxSection[0] ; }
	  G4double GetMeanCerenkovLoss() const {return fIntegralCerenkov[0] ; }
	  G4double GetMeanPlasmonLoss() const {return fIntegralPlasmon[0] ; }

	  G4double GetNormalizationCof() const { return fNormalizationCof ; }
          
	  inline G4double GetPAItable(G4int i,G4int j) const ;

          inline G4double    GetLorentzFactor(G4int i) const ;
	  	  
	  inline G4double GetSplineEnergy(G4int i) const ;
	  
	  inline G4double GetIntegralPAIxSection(G4int i) const ;
	  inline G4double GetIntegralCerenkov(G4int i) const ;
	  inline G4double GetIntegralPlasmon(G4int i) const ;

protected :

private :

// Local class constants
 
static const G4double fDelta ; // energy shift from interval border = 0.001
static const G4double fError ; // error in lin-log approximation = 0.005

static       G4int fNumberOfGammas ;         // = 111 ;
static const G4double fLorentzFactor[112] ;  //  static gamma array

static 
const G4int fRefGammaNumber  ; // The number of gamma for creation of spline (15)

G4int    fIntervalNumber  ;    //  The number of energy intervals
G4double fNormalizationCof ;   // Normalization cof for PhotoAbsorptionXsection
// G4double fBetaGammaSq ;        // (beta*gamma)^2

G4double fDensity ;            // Current density
G4double fElectronDensity ;    // Current electron (number) density
G4int    fSplineNumber ;       // Current size of spline

// Arrays of Sandia coefficients

G4double* fEnergyInterval ;
G4double* fA1 ; 
G4double* fA2 ;
G4double* fA3 ; 
G4double* fA4 ;

static
const G4int   fMaxSplineSize  ;          // Max size of output splain arrays = 500
/* ******************
G4double*          fSplineEnergy ;   // energy points of splain
G4double* fRePartDielectricConst ;   // Real part of dielectric const
G4double* fImPartDielectricConst ;   // Imaginary part of dielectric const
G4double*          fIntegralTerm ;   // Integral term in PAI cross section
G4double*        fDifPAIxSection ;   // Differential PAI cross section
G4double*   fIntegralPAIxSection ;   // Integral PAI cross section  ?
*/ ///////////////
G4double          fSplineEnergy[500] ;   // energy points of splain
G4double fRePartDielectricConst[500] ;   // Real part of dielectric const
G4double fImPartDielectricConst[500] ;   // Imaginary part of dielectric const
G4double          fIntegralTerm[500] ;   // Integral term in PAI cross section
G4double        fDifPAIxSection[500] ;   // Differential PAI cross section
G4double          fdNdxCerenkov[500] ;   // dNdx of Cerenkov collisions
G4double          fdNdxPlasmon[500] ;   // dNdx of Plasmon collisions

G4double   fIntegralPAIxSection[500] ;   // Integral PAI cross section  ?
G4double   fIntegralCerenkov[500] ;   // Integral Cerenkov N>omega  ?
G4double   fIntegralPlasmon[500] ;   // Integral Plasmon N>omega  ?

G4double fPAItable[500][112] ; // Output array

} ;    

////////////////  Inline methods //////////////////////////////////
//


inline G4double G4PAIxSection::GetPAItable(G4int i, G4int j) const
{
   return fPAItable[i][j] ;
}

inline G4double G4PAIxSection::GetLorentzFactor(G4int j) const
{
   return fLorentzFactor[j] ;
}

inline G4double G4PAIxSection::GetSplineEnergy(G4int i) const 
{
   if(i < 1 || i > fSplineNumber)
   {
      G4Exception("Invalid argument in G4PAIxSection::GetSplineEnergy");
   }
   return fSplineEnergy[i] ;
}
	  
inline G4double G4PAIxSection::GetIntegralPAIxSection(G4int i) const 
{
   if(i < 1 || i > fSplineNumber)
   {
    G4Exception("Invalid argument in G4PAIxSection::GetIntegralPAIxSection");
   }
   return fIntegralPAIxSection[i] ;
}

inline G4double G4PAIxSection::GetIntegralCerenkov(G4int i) const 
{
   if(i < 1 || i > fSplineNumber)
   {
    G4Exception("Invalid argument in G4PAIxSection::GetIntegralCerenkov");
   }
   return fIntegralCerenkov[i] ;
}

inline G4double G4PAIxSection::GetIntegralPlasmon(G4int i) const 
{
   if(i < 1 || i > fSplineNumber)
   {
    G4Exception("Invalid argument in G4PAIxSection::GetIntegralPlasmon");
   }
   return fIntegralPlasmon[i] ;
}

#endif   

// -----------------   end of G4PAIxSection header file    -------------------
