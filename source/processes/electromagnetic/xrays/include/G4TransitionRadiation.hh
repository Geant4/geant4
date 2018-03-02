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
// $Id: G4TransitionRadiation.hh 108508 2018-02-15 15:54:35Z gcosmo $
//
// G4TransitionRadiation  -- header file
//
// Class for description of  transition radiation generated
// by  charged particle crossed interface between material 1
// and material 2 (1 -> 2). Transition radiation could be of kind:
// - optical back
// - optical forward
// - X-ray   forward (for relativistic case Tkin/mass >= 10^2)
//
// GEANT 4 class header file --- Copyright CERN 1995
// CERB Geneva Switzerland
//
// for information related to this code, please, contact
// CERN, CN Division, ASD Group
// History:
// 18.12.97, V. Grichine (Vladimir.Grichine@cern.ch)
// 02.02.00, V.Grichine, new data fEnergy and fVarAngle for double
//                       numerical integration in inherited classes
// 03.06.03, V.Ivanchenko fix compilation warnings
// 28.07.05, P.Gumplinger add G4ProcessType to constructor

#ifndef G4TransitionRadiation_h
#define G4TransitionRadiation_h


#include "G4VDiscreteProcess.hh"
#include "G4Material.hh"

class G4TransitionRadiation : public   G4VDiscreteProcess
{
public:

  explicit G4TransitionRadiation( const G4String& processName = "TR",
			 G4ProcessType type = fElectromagnetic) ;

  virtual ~G4TransitionRadiation() ;

  // Methods

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  virtual G4double GetMeanFreePath(const G4Track&, G4double,
			   G4ForceCondition* condition) override;

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, 
                                          const G4Step&) override;

  virtual
  G4double SpectralAngleTRdensity( G4double energy,
                                 G4double varAngle ) const = 0 ;

  G4double IntegralOverEnergy( G4double energy1,
			       G4double energy2,
			       G4double varAngle     ) const ;

  G4double IntegralOverAngle( G4double energy,
			      G4double varAngle1,
			      G4double varAngle2     ) const ;

  G4double AngleIntegralDistribution( G4double varAngle1,
				      G4double varAngle2     ) const ;

  G4double EnergyIntegralDistribution( G4double energy1,
				       G4double energy2     )  const   ;



  // Access functions

protected :

  G4int fMatIndex1 ;                   // index of the 1st material
  G4int fMatIndex2 ;                   // index of the 2nd material

  // private :

  G4double fGamma ;
  G4double fEnergy ;
  G4double fVarAngle ;

  // Local constants
  static const G4int fSympsonNumber ; // Accuracy of Sympson integration 10
  static const G4int fGammaNumber   ; // = 15
  static const G4int fPointNumber   ; // = 100

  G4double fMinEnergy ;                //  min TR energy
  G4double fMaxEnergy ;                //  max TR energy
  G4double fMaxTheta  ;                //  max theta of TR quanta

  G4double fSigma1 ;                   // plasma energy Sq of matter1
  G4double fSigma2 ;                   // plasma energy Sq of matter2

private:

// Operators
  G4TransitionRadiation(const G4TransitionRadiation& right) = delete;
  G4TransitionRadiation& 
    operator=(const G4TransitionRadiation& right) = delete;

};

#endif   // G4TransitionRadiation_h
