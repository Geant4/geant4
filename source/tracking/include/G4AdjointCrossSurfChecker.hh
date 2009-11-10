/////////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointCrossSurfChecker
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 January 2007 creation by L. Desorgher 
//		-01/11/2009 Some cleaning and adding of documentation  for the first Release in the Geant4 toolkit, L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		This class is responsible for checking the crossing of a surface that could be the  external boundary of a volume  or the external surface of a sphere.
//		It is used to check if an adjoint particle reaches the external surface or reenter the sensitive region delimited by the adjoint source.   
//

#ifndef G4AdjointCrossSurfChecker_h
#define G4AdjointCrossSurfChecker_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4Step;


class G4AdjointCrossSurfChecker 
{ 
 public:
  
   G4AdjointCrossSurfChecker();
   virtual ~G4AdjointCrossSurfChecker();
   
   static G4AdjointCrossSurfChecker* GetInstance();
  
 public:
   
   bool CrossingASphere(const G4Step* aStep,G4double sphere_radius, G4ThreeVector sphere_center,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   bool GoingInOrOutOfaVolume(const G4Step* aStep,G4String volume_name, G4double& cos_to_surface, G4bool& GoingIn);
   bool GoingInOrOutOfaVolumeByExternalSurface(const G4Step* aStep,G4String volume_name,G4String mother_log_vol_name, G4double& cos_to_surface, G4bool& GoingIn);
   
 
   bool CrossingAGivenRegisteredSurface(const G4Step* aStep,G4String surface_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   bool CrossingAGivenRegisteredSurface(const G4Step* aStep, int ind,G4ThreeVector& crossing_pos,   G4double& cos_to_surface, G4bool& GoingIn);
   bool CrossingOneOfTheRegisteredSurface(const G4Step* aStep,G4String& surface_name,G4ThreeVector& crossing_pos,G4double& cos_to_surface, G4bool& GoingIn);
   bool CrossingAnInterfaceBetweenTwoVolumes(const G4Step* aStep,G4String vol1_name,G4String vol2_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   
   bool AddaSphericalSurface(G4String SurfaceName, G4double radius, G4ThreeVector pos,G4double& area);
   bool AddaSphericalSurfaceWithCenterAtTheCenterOfAVolume(G4String SurfaceName, G4double radius, G4String  volume_name, G4ThreeVector& center,  G4double& area);
   bool AddanExternalSurfaceOfAvolume(G4String SurfaceName, G4String volume_name,G4double& area);
   bool AddanInterfaceBetweenTwoVolumes(G4String SurfaceName, G4String volume_name1, G4String volume_name2,G4double& area);
   void ClearListOfSelectedSurface();
    
  
   
   		       
 private: 
   
   int FindRegisteredSurface(G4String name);
   
 private: 
   static G4AdjointCrossSurfChecker* instance;
  
   std::vector<G4String> ListOfSurfaceName;
   std::vector<G4String> ListOfSurfaceType;
   std::vector<G4double> ListOfSphereRadius;
   std::vector<G4ThreeVector> ListOfSphereCenter;
   std::vector<G4String> ListOfVol1Name;
   std::vector<G4String> ListOfVol2Name;
   std::vector<G4double> AreaOfSurface;

};

#endif
