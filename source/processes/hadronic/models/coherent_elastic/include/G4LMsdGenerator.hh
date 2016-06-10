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
// $Id$
//
// Author: Vladimir Grichine, ~25 October 2014 
//
// Class Description
// Final state production model for theoretical models of hadron inelastic
// scattering in geant4;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Note: This class is part of an implementation framework. You need to
// register corresponding high energy generators and transport codes to 
// fill it with life; decay of strong resonances is done directly,
// in case there is no residual nucleus. 
// Class Description - End
//
// History:
//
// ~25.10.14 V. Grichine: p and kinematics tests
// 20.11.14 V. Grichine: n, pi+-, K+- added
//
//



#ifndef G4LMsdGenerator_h
#define G4LMsdGenerator_h 1

#include "G4HadronicInteraction.hh"
#include "G4HadFinalState.hh"

// class G4ParticleDefinition;

class G4LMsdGenerator : public G4HadronicInteraction 

{
  public:

  G4LMsdGenerator(const G4String& name = "LMsdGenerator");
  ~G4LMsdGenerator();

  private:

  G4LMsdGenerator(const G4LMsdGenerator &right);
  const G4LMsdGenerator & operator=(const G4LMsdGenerator &right);
  int operator == (const G4LMsdGenerator &right) const;
  int operator != (const G4LMsdGenerator &right) const;

  public:

  G4bool IsApplicable(const G4HadProjectile & thePrimary, 
                                        G4Nucleus & theNucleus);

  G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, 
                                        G4Nucleus & theNucleus);

  G4double SampleMx(const G4HadProjectile* aParticle );

  G4double SampleT( const G4HadProjectile* aParticle, G4double Mx );

  void ModelDescription(std::ostream& outFile) const;

  private: 

  // G4ParticleDefinition* fParticle;

  G4int fPDGencoding;

  static const G4double fMxBdata[23][2];
  static const G4double  fProbMx[60][2];
};


#endif


