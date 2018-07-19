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
// $Id: G4GammaTransition.hh 85659 2014-11-03 10:59:10Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4GammaTransition
//
//      Author:        Vladimir Ivanchenko
// 
//      Creation date: 27 February 2015
//
//      Modifications: 
// 05.11.2016 J. Detviller added sampling of corelation between
//            spin of a nucleus and gamma direction
//
//      
// -------------------------------------------------------------------

#ifndef G4GAMMATRANSITION_HH
#define G4GAMMATRANSITION_HH 1

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4PolarizationTransition.hh"

class G4GammaTransition 
{
public:

  explicit G4GammaTransition();

  virtual ~G4GammaTransition();
  
  virtual G4Fragment* SampleTransition(G4Fragment* nucleus,
				       G4double newExcEnergy,
                                       G4double mpRatio,
                                       G4int  JP1,
                                       G4int  JP2,
                                       G4int  MP,
                                       G4int  shell,
                                       G4bool isDiscrete,
                                       G4bool isGamma);

  virtual void SampleDirection(G4Fragment* nuc, G4double ratio, 
			       G4int twoJ1, G4int twoJ2, G4int mp);

  inline void SetPolarizationFlag(G4bool val) { polarFlag = val; };

  inline void SetTwoJMAX(G4int val) { fTwoJMAX = val; };

  inline void SetVerbose(G4int val) { fVerbose = val; fPolTrans.SetVerbose(val); };

private:  

  G4GammaTransition(const G4GammaTransition &right) = delete;
  const G4GammaTransition& operator=(const G4GammaTransition &right) = delete;
  G4bool operator==(const G4GammaTransition &right) const = delete;
  G4bool operator!=(const G4GammaTransition &right) const = delete;
 
  G4bool polarFlag;

protected:

  G4ThreeVector fDirection;
  G4PolarizationTransition fPolTrans;
  G4int fTwoJMAX;
  G4int fVerbose;
};


#endif
