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
// $Id: G4InitXscPAI.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// 
// G4InitXscPAI.hh -- header file
//
// History:
//
// 02.04.04, V. Grichine: 1st version based on G4PAIxSection class

#ifndef G4INITXSCPAI_HH
#define G4INITXSCPAI_HH

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4OrderedTable.hh"
#include "G4PhysicsLogVector.hh"

class G4MaterialCutsCouple;
class G4SandiaTable;

class G4InitXscPAI
{
public:
	  // Constructors
  explicit G4InitXscPAI( const G4MaterialCutsCouple* matCC);

  virtual ~G4InitXscPAI() ;

  // Methods
  // General control functions

  void KillCloseIntervals();

  void Normalisation();


  // Physical methods


  G4double RutherfordIntegral( G4int intervalNumber,
	                               G4double limitLow,
				       G4double limitHigh     ) ;

  G4double IntegralTerm(G4double omega);

  G4double ImPartDielectricConst( G4int intervalNumber,
	                                  G4double energy        ) ;

  G4double RePartDielectricConst(G4double energy) ;

  G4double ModuleSqDielectricConst( G4int intervalNumber,
	                                  G4double energy        ) ;

  G4double DifPAIxSection( G4double omega ) ;
  G4double DifPAIdEdx( G4double omega ) ;

  G4double PAIdNdxCherenkov( G4double omega ) ;

  G4double PAIdNdxPlasmon( G4double omega ) ;

  void     IntegralPAIxSection(G4double bg2, G4double Tmax) ;
  void     IntegralCherenkov(G4double bg2, G4double Tmax) ;
  void     IntegralPlasmon(G4double bg2, G4double Tmax) ;

  void      IntegralPAIdEdx(G4double bg2, G4double Tmax) ;


  G4double GetPhotonLambda( G4double omega ) ;


  G4double GetStepEnergyLoss( G4double step ) ;
  G4double GetStepCerenkovLoss( G4double step ) ;
  G4double GetStepPlasmonLoss( G4double step ) ;

  // Inline access functions


  G4int GetIntervalNumber() const { return fIntervalNumber ; }
  G4int GetBinPAI() const { return fPAIbin ; }

  G4double GetNormalizationCof() const { return fNormalizationCof ; }

  G4double GetMatSandiaMatrix(G4int i, G4int j) const
          { return (*(*fMatSandiaMatrix)[i])[j]; }

  G4PhysicsLogVector* GetPAIxscVector() const { return fPAIxscVector;}
  G4PhysicsLogVector* GetPAIdEdxVector() const { return fPAIdEdxVector;}
  G4PhysicsLogVector* GetPAIphotonVector() const { return fPAIphotonVector;}
  G4PhysicsLogVector* GetPAIelectronVector() const { return fPAIelectronVector;}
  G4PhysicsLogVector* GetChCosSqVector() const { return fChCosSqVector;}
  G4PhysicsLogVector* GetChWidthVector() const { return fChWidthVector;}

protected :

private :

  G4InitXscPAI & operator=(const G4InitXscPAI &right) = delete;
  G4InitXscPAI(const G4InitXscPAI&) = delete;

  // Local class constants

  static const G4double fDelta ; // energy shift from interval border = 0.001
  static const G4int fPAIbin;
  static const G4double fSolidDensity; // ~the border between gases and solids

  G4int    fIntervalNumber;    //  The number of energy intervals
  G4double fNormalizationCof;   // Normalization cof for PhotoAbsorptionXsection
  G4int    fCurrentInterval;
  G4int    fIntervalTmax;
  G4double fBetaGammaSq ;        // (beta*gamma)^2
  G4double fTmax;
  G4double fDensity ;            // Current density
  G4double fElectronDensity ;    // Current electron (number) density

  // Arrays of Sandia coefficients

  G4OrderedTable* fMatSandiaMatrix;
  G4SandiaTable*  fSandia;

  // vectors of integral cross-sections
  
  G4PhysicsLogVector* fPAIxscVector;
  G4PhysicsLogVector* fPAIdEdxVector;
  G4PhysicsLogVector* fPAIphotonVector;
  G4PhysicsLogVector* fPAIelectronVector;
  G4PhysicsLogVector* fChCosSqVector;
  G4PhysicsLogVector* fChWidthVector;
  
};    

#endif   

//
//
/////////////////   end of G4InitXscPAI header file    //////////////////////
