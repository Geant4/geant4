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
#ifndef G4ParticleHPJENDLHEData_h
#define G4ParticleHPJENDLHEData_h 1

// Class Description
// Cross-section data set for a high precision (based on JENDL_HE evaluated data
// libraries) description of elastic scattering 20 MeV ~ 3 GeV; 
// Class Description - End

// 15-Nov-06 First Implementation is done by T. Koi (SLAC/SCCS)
// 080422 Add IsZAApplicable method (return false) by T. Koi
//

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include <map> 

class G4ParticleHPJENDLHEData : public G4VCrossSectionDataSet
{
   public:
   
   G4ParticleHPJENDLHEData();
   G4ParticleHPJENDLHEData( G4String , G4ParticleDefinition* );

   
   ~G4ParticleHPJENDLHEData();
   
   G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

   G4bool IsZAApplicable( const G4DynamicParticle* , G4double /*ZZ*/, G4double /*AA*/)
   { return false; }

   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);

   void BuildPhysicsTable(const G4ParticleDefinition&);

   void DumpPhysicsTable(const G4ParticleDefinition&);
   
   private:
   
      std::vector< G4bool > vElement;

      std::map< G4int , std::map< G4int , G4PhysicsVector* >* > mIsotope; 

      G4bool isThisInMap ( G4int , G4int );
      G4bool isThisNewIsotope ( G4int z , G4int a ) { return !( isThisInMap( z , a ) ); };
      G4PhysicsVector* readAFile ( std::fstream* );
      void registAPhysicsVector ( G4int , G4int , G4PhysicsVector* );

      G4double getXSfromThisIsotope ( G4int , G4int , G4double );

      G4String reactionName;
      G4String particleName;
};

#endif
