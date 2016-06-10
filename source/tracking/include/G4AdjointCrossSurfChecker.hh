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
// $Id: G4AdjointCrossSurfChecker.hh 66872 2013-01-15 01:25:57Z japost $
//
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
//		- 01/11/2009 Some cleaning and adding of documentation for the
//                           first Release in the Geant4 toolkit, L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		This class is responsible for checking the crossing of a surface
//              that could be the  external boundary of a volume  or the external
//              surface of a sphere.
//		It is used to check if an adjoint particle reaches the external
//              surface or reenter the sensitive region delimited by the adjoint
//              source.   
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
  
   static G4AdjointCrossSurfChecker* GetInstance();
  
 public:
   
   G4bool CrossingASphere(const G4Step* aStep,G4double sphere_radius, G4ThreeVector sphere_center,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   G4bool GoingInOrOutOfaVolume(const G4Step* aStep,const G4String& volume_name, G4double& cos_to_surface, G4bool& GoingIn);
   G4bool GoingInOrOutOfaVolumeByExtSurface(const G4Step* aStep,const G4String& volume_name,const G4String& mother_log_vol_name, G4double& cos_to_surface, G4bool& GoingIn);

   G4bool CrossingAGivenRegisteredSurface(const G4Step* aStep,const G4String& surface_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   G4bool CrossingAGivenRegisteredSurface(const G4Step* aStep, int ind,G4ThreeVector& crossing_pos,   G4double& cos_to_surface, G4bool& GoingIn);
   G4bool CrossingOneOfTheRegisteredSurface(const G4Step* aStep,G4String& surface_name,G4ThreeVector& crossing_pos,G4double& cos_to_surface, G4bool& GoingIn);
   G4bool CrossingAnInterfaceBetweenTwoVolumes(const G4Step* aStep,const G4String& vol1_name,const G4String& vol2_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn);
   
   G4bool AddaSphericalSurface(const G4String& SurfaceName, G4double radius, G4ThreeVector pos,G4double& area);
   G4bool AddaSphericalSurfaceWithCenterAtTheCenterOfAVolume(const G4String& SurfaceName, G4double radius, const G4String& volume_name, G4ThreeVector& center,  G4double& area);
   G4bool AddanExtSurfaceOfAvolume(const G4String& SurfaceName,const G4String& volume_name,G4double& area);
   G4bool AddanInterfaceBetweenTwoVolumes(const G4String& SurfaceName, const G4String& volume_name1, const G4String& volume_name2,G4double& area);
   void ClearListOfSelectedSurface();
   		       
 private: 
   
   G4AdjointCrossSurfChecker();
  ~G4AdjointCrossSurfChecker();
   
   G4int FindRegisteredSurface(const G4String& name);
   
 private: 
   static G4ThreadLocal G4AdjointCrossSurfChecker* instance;
  
   std::vector<G4String> ListOfSurfaceName;
   std::vector<G4String> ListOfSurfaceType;
   std::vector<G4double> ListOfSphereRadius;
   std::vector<G4ThreeVector> ListOfSphereCenter;
   std::vector<G4String> ListOfVol1Name;
   std::vector<G4String> ListOfVol2Name;
   std::vector<G4double> AreaOfSurface;

};

#endif

