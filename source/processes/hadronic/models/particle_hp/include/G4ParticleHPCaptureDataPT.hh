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
// -------------------------------------------------------------------
//
//      Geant4 header file 
//
//      File name: G4ParticleHPCaptureDataPT.hh
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//	         Loic Thulliez (CEA France)
// 
//      Creation date: 4 June 2024
//
//      Description: Class for utilization of cross-sections from
//                   probability tables in the unresolved resonance region
//                   for capture channel.
//                   Cross-section data set for a high precision
//                   (based on evaluated data libraries) description of
//                   neutron Capture scattering below 20 MeV.
//                   To be used in your physics list in case you need
//                   this physics.
//                   In this case you want to register an object of this
//                   class with the corresponding process.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//
#ifndef G4ParticleHPCaptureDataPT_h
#define G4ParticleHPCaptureDataPT_h 1

#include "G4VCrossSectionDataSet.hh"
#include <vector>

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;


class G4ParticleHPCaptureDataPT : public G4VCrossSectionDataSet {
  public:
    G4ParticleHPCaptureDataPT();
    ~G4ParticleHPCaptureDataPT();

    void BuildPhysicsTable( const G4ParticleDefinition& );
    G4bool IsIsoApplicable( const G4DynamicParticle* , G4int /*Z*/ , G4int /*A*/ ,
                            const G4Element* /*elm*/ , const G4Material* /*mat*/ );
    G4double GetIsoCrossSection( const G4DynamicParticle* , G4int /*Z*/ , G4int /*A*/ ,
                                 const G4Isotope* /*iso*/ , const G4Element* /*elm*/ , const G4Material* /*mat*/ );

    void SetVerboseLevel( G4int );
    G4int GetVerboseLevel() const;
    virtual void CrossSectionDescription( std::ostream& ) const;

   private:
     std::vector< std::pair< G4double, G4double > >* URRlimits;
};

#endif
