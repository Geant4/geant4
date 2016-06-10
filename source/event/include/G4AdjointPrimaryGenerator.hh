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
// $Id: G4AdjointPrimaryGenerator.hh 86965 2014-11-21 11:48:22Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointPrimaryGenerator
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//		November 2009 creation by L. Desorgher, Splitting of G4AdjointPrimaryGeneratorAction in two classes  G4AdjointPrimaryGeneratorAction and G4AdjointPrimaryGenerator  		
//
//-------------------------------------------------------------
//	Documentation:
//		This class represents the Primary Generator that generate vertex (energy,position and direction) of primary adjoint particles.
//		It is used by  G4AdjointPrimaryGeneratorAction. If the adjoint source is selected by the user as being on the external boundary of a volume
//		it uses the class G4AdjointPosOnPhysVolGenerator to generate the vertex positions and directions. Otherwise  G4SingleParticleSource is used. 
//		
//		
//
#ifndef G4AdjointPrimaryGenerator_h
#define G4AdjointPrimaryGenerator_h 1
#include "globals.hh"
#include"G4ThreeVector.hh"
#include <vector>
#include <map>
#include <iterator>
#include"G4PhysicsOrderedFreeVector.hh"


class G4AdjointPosOnPhysVolGenerator;
class G4Event;
class G4SingleParticleSource;
class G4ParticleDefinition;
class G4Navigator;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class G4AdjointPrimaryGenerator
{ public:
    G4AdjointPrimaryGenerator();    
   ~G4AdjointPrimaryGenerator();

  public: //public methods
    
    void GenerateAdjointPrimaryVertex(G4Event* anEvt,G4ParticleDefinition* adj_part,G4double E1,G4double E2);
    void GenerateFwdPrimaryVertex(G4Event* anEvt,G4ParticleDefinition* adj_part,G4double E1,G4double E2);
    void SetSphericalAdjointPrimarySource(G4double radius, G4ThreeVector pos);
    void SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(const G4String& volume_name);
    void ComputeAccumulatedDepthVectorAlongBackRay(G4ThreeVector glob_pos,
                                                   G4ThreeVector direction,
                                                   G4double ekin,
                                                   G4ParticleDefinition* aPartDef);
    G4double SampleDistanceAlongBackRayAndComputeWeightCorrection(G4double& weight_corr);


  private: //attributes

    //The class responsible for the random generation of  positions and direction of primaries for adjoint source set on the external surface of
    //a G4 volume
    G4AdjointPosOnPhysVolGenerator* theG4AdjointPosOnPhysVolGenerator;
    
    G4SingleParticleSource* theSingleParticleSource;
    
    //Type of adjoint source
    //--------------------
    G4String type_of_adjoint_source; //Spherical ExtSurfaceOfAVolume
    G4double radius_spherical_source;
    G4ThreeVector center_spherical_source;
    G4Navigator* fLinearNavigator;
    G4PhysicsOrderedFreeVector* theAccumulatedDepthVector;
    //G4PhysicsOrderedFreeVector* theAccumulatedCSDepthVector;

    //Disable copy constructor and assignement operator
    G4AdjointPrimaryGenerator(const G4AdjointPrimaryGenerator&);
    G4AdjointPrimaryGenerator& operator=(const G4AdjointPrimaryGenerator&);



};
#endif

