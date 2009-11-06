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

class G4AdjointPosOnPhysVolGenerator;
class G4Event;
class G4SingleParticleSource;
class G4ParticleDefinition;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class G4AdjointPrimaryGenerator 
{
  public: //constructor, destructor
    G4AdjointPrimaryGenerator();    
   ~G4AdjointPrimaryGenerator();

  public: //public methods
    
    void GenerateAdjointPrimaryVertex(G4Event* anEvt,G4ParticleDefinition* adj_part,G4double E1,G4double E2);
    void SetSphericalAdjointPrimarySource(G4double radius, G4ThreeVector pos);
    void SetAdjointPrimarySourceOnAnExternalSurfaceOfAVolume(G4String volume_name);
    
  
  private: //attributes
    
    
    
    //The class responsible for the random generation of  positions and direction of primaries for adjoint source set on the external surface of
    //a G4 volume
    G4AdjointPosOnPhysVolGenerator* theG4AdjointPosOnPhysVolGenerator;
    
    
    G4SingleParticleSource* theSingleParticleSource;
    
    //Type of adjoint source
    //--------------------
    G4String type_of_adjoint_source; //Spherical ExternalSurfaceOfAVolume
    G4double radius_spherical_source;
    G4ThreeVector center_spherical_source;
  
   
};
#endif


