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
#ifndef G4LENDCombinedModel_h
#define G4LENDCombinedModel_h 1

// Class Description
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 170912 First implementation done by T. Koi (SLAC/EPP)

#include "G4LENDModel.hh"

class G4LENDCombinedCrossSection;

class G4LENDElastic;
class G4LENDInelastic;
class G4LENDCapture;
class G4LENDFission;

class G4LENDCombinedModel : public G4LENDModel
{
   public: 
      G4LENDCombinedModel( G4ParticleDefinition* pd );
      ~G4LENDCombinedModel(){;};

      void BuildPhysicsTable(const G4ParticleDefinition&);

      G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );

      G4bool HasData( const G4DynamicParticle* , G4int iZ , G4int iA , G4int iM, 
                      const G4Isotope* , const G4Element* , const G4Material* );
     
   private: 
      G4LENDCombinedCrossSection* crossSection;
      G4LENDElastic* elastic;
      G4LENDInelastic* inelastic;
      G4LENDCapture* capture;
      G4LENDFission* fission;
      G4LENDModel* channels[4];
};

#endif
