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
 
#ifndef G4LENDFission_h
#define G4LENDFission_h 1

// Class Description
// Manager of LEND (Low Energy Nuclear Data) target (nucleus) 
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDModel.hh"

class G4LENDFission : public G4LENDModel
{

   public: 
  
     G4LENDFission( G4ParticleDefinition* pd )
     :G4LENDModel( "LENDFission" ) 
     { 
        proj = pd; 
       
//        theModelName = "LENDFission for "; 
//        theModelName += proj->GetParticleName(); 
        create_used_target_map();
     };
  
     ~G4LENDFission(){;};
     virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;
     G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );
};

#endif
