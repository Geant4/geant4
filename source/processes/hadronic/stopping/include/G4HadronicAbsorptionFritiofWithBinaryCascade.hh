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
//---------------------------------------------------------------------------
//
// ClassName:  G4HadronicAbsorptionFritiofWithBinaryCascade
//
// Author:     Alberto Ribon
//
// Date:       July 2019
//
// Modified:
//
// Class Description:
//
// Intermediate base class for hadronic absorption at rest using
// Fritiof (FTF) coupled with Binary Cascade (BIC) for the nuclear
// de-excitation. 
//
// Note: it is applicable only for anti_proton and anti_neutron,
//       but not for light anti-ion because the propagate method
//       for nucleus-nucleus interaction:
//         G4VIntraNuclearTransportModel::Propagate()
//       is not implemented.
//
//---------------------------------------------------------------------------

#ifndef G4HadronicAbsorptionFritiofWithBinaryCascade_h
#define G4HadronicAbsorptionFritiofWithBinaryCascade_h 1

#include "globals.hh"
#include "G4HadronStoppingProcess.hh"

class G4ParticleDefinition;

class G4HadronicAbsorptionFritiofWithBinaryCascade : public G4HadronStoppingProcess { 

public:
  G4HadronicAbsorptionFritiofWithBinaryCascade( G4ParticleDefinition* pdef = 0 ); 
  virtual ~G4HadronicAbsorptionFritiofWithBinaryCascade();
  
  G4bool IsApplicable( const G4ParticleDefinition& );

  void ProcessDescription( std::ostream& outFile ) const;

private:
  // hide assignment operator as private 
  G4HadronicAbsorptionFritiofWithBinaryCascade& operator=( const G4HadronicAbsorptionFritiofWithBinaryCascade& );
  G4HadronicAbsorptionFritiofWithBinaryCascade( const G4HadronicAbsorptionFritiofWithBinaryCascade& );
  
  G4ParticleDefinition* pdefApplicable;
};

#endif

