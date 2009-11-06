/////////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPosOnPhysVolGenerator
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1st June 2006 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		This class is responsible for the generation of primary adjoint particle on the external surface of a user selected volume.
//		The particle are generated uniformly on the surface with the angular distribution  set to a cosine law relative to normal of the surface.
//		It is equivalent to  the flux going in  from the surface if an  isotropic flux is considered outside. 
//		It uses ray tracking technique and can be applied to all kind of convex volume. Uisng the ray tracking technique the area 
//		of the external surface is also computed. The area is needed to fix the weight of the primary adjoint particle.  
//		At the time of the development of this class, generation of particle on volume surface and computation of surface was limited in G4, 
//		therfore the general ray tracking technique was adopted. It could be now (2009) that direct method of G4VSolid could be used instead. To be checked! 
//
//		
//
#ifndef G4AdjointPosOnPhysVolGenerator_h
#define G4AdjointPosOnPhysVolGenerator_h 1



class G4VSolid;
#include "G4VPhysicalVolume.hh"
#include "G4AffineTransform.hh"

#include "G4ThreeVector.hh"
class G4AdjointPosOnPhysVolGenerator 
///////////////////////
{

//--------
  public: //without description
//--------

// Constructor/Destructor
   
   static  G4AdjointPosOnPhysVolGenerator* GetInstance();
   
   ~G4AdjointPosOnPhysVolGenerator();
   
//--------
  public:  //public methods
//--------
  G4VPhysicalVolume* DefinePhysicalVolume(G4String aName);
  void DefinePhysicalVolume1(G4String aName);
  double ComputeAreaOfExternalSurface();
  double ComputeAreaOfExternalSurface(G4int NStat);
  double ComputeAreaOfExternalSurface(double epsilon);
  double ComputeAreaOfExternalSurface(G4VSolid* aSolid);
  double ComputeAreaOfExternalSurface(G4VSolid* aSolid,G4int NStat);
  double ComputeAreaOfExternalSurface(G4VSolid* aSolid,double epsilon);
 
  void GenerateAPositionOnTheExternalSurfaceOfASolid(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction);
  void GenerateAPositionOnTheExternalSurfaceOfTheSolid(G4ThreeVector& p, G4ThreeVector&  direction);
  void GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(G4ThreeVector& p, G4ThreeVector&  direction);
  void GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(G4ThreeVector& p, G4ThreeVector&  direction,
  										G4double& costh_to_normal);

  //inline public methods
   
  inline void SetSolid(G4VSolid* aSolid){theSolid=aSolid;}
  inline G4double GetAreaOfExternalSurfaceOfThePhysicalVolume(){return AreaOfExternalSurfaceOfThePhysicalVolume;}
  inline G4double GetCosThDirComparedToNormal(){return CosThDirComparedToNormal;}
  
//---------   
   private:   //private methods
//---------  
   G4AdjointPosOnPhysVolGenerator();
   double ComputeAreaOfExternalSurfaceStartingFromSphere(G4VSolid* aSolid,G4int NStat);
   double ComputeAreaOfExternalSurfaceStartingFromBox(G4VSolid* aSolid,G4int NStat);
   void GenerateAPositionOnASolidBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction);
   double GenerateAPositionOnASphereBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction);
   double GenerateAPositionOnABoxBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction);
   void ComputeTransformationFromPhysVolToWorld();

//---------   
   private: //attributes
//---------   
   static G4AdjointPosOnPhysVolGenerator* theInstance;
   G4VSolid* theSolid;
   G4VPhysicalVolume* thePhysicalVolume;
   G4int NStat;
   double epsilon;
   bool UseSphere;
   G4String ModelOfSurfaceSource;
   double ExternalSourceRadius;
   double ExternalSourceDx,ExternalSourceDy,ExternalSourceDz ;
   G4AffineTransform theTransformationFromPhysVolToWorld;
   double AreaOfExternalSurfaceOfThePhysicalVolume;
   double CosThDirComparedToNormal;
   

 
   	

};



#endif
