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
// $Id: G4AdjointSimManager.cc,v 1.2 2009-11-18 18:02:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointCrossSurfChecker
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////

#include "G4AdjointSimManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include "G4UserEventAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserRunAction.hh"

#include "G4AdjointPrimaryGeneratorAction.hh"
#include "G4AdjointSteppingAction.hh"
#include "G4AdjointStackingAction.hh"

#include "G4AdjointSimMessenger.hh"

#include "G4AdjointCrossSurfChecker.hh"

#include "G4ParticleTable.hh"
#include "G4PhysicsLogVector.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSimManager* G4AdjointSimManager::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSimManager::G4AdjointSimManager()
{ 
 //Create adjoint actions;
 //----------------------
 
 theAdjointRunAction = 0; 
 theAdjointPrimaryGeneratorAction = new G4AdjointPrimaryGeneratorAction(); 
 theAdjointSteppingAction = new G4AdjointSteppingAction();
 theAdjointEventAction = 0;
 theAdjointTrackingAction = 0;	
 theAdjointStackingAction = new G4AdjointStackingAction();

 //Create messenger
 //----------------
  theMessenger = new G4AdjointSimMessenger(this);
  
  user_action_already_defined=false;
  use_user_StackingAction = false; 

  fUserTrackingAction= 0;
  fUserEventAction= 0; 
  fUserSteppingAction= 0;
  fUserPrimaryGeneratorAction= 0;
  fUserRunAction= 0;
  fUserStackingAction= 0;
  
  adjoint_sim_mode = false;
  
  normalisation_mode=3;
  
  nb_nuc=1.;
  
  welcome_message =true;
  
  /*electron_last_weight_vector = new G4PhysicsLogVector(1.e-20,1.e20,400);
  proton_last_weight_vector  = new  G4PhysicsLogVector(1.e-20,1.e20,400);
  gamma_last_weight_vector = new    G4PhysicsLogVector(1.e-20,1.e20,400);*/ 
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSimManager::~G4AdjointSimManager()
{ 
  if (theAdjointRunAction) delete theAdjointRunAction;
  if (theAdjointPrimaryGeneratorAction) delete theAdjointPrimaryGeneratorAction;
  if (theAdjointSteppingAction) delete theAdjointSteppingAction;
  if (theAdjointEventAction) delete theAdjointEventAction;
  if (theAdjointTrackingAction) delete theAdjointTrackingAction;
  if (theAdjointStackingAction) delete theAdjointStackingAction;
  if (theMessenger) delete theMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSimManager* G4AdjointSimManager::GetInstance()
{
  if (instance == 0) instance = new G4AdjointSimManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::RunAdjointSimulation(G4int nb_evt)
{ 
  if (welcome_message) {
  	G4cout<<"****************************************************************"<<std::endl;
	G4cout<<"*** Geant4 Reverse/Adjoint Monte Carlo mode		      ***"<<std::endl;
	G4cout<<"*** Author:	L.Desorgher				      ***"<<std::endl;
	G4cout<<"*** Company:	SpaceIT GmbH, Bern, Switzerland 	      ***"<<std::endl;
	G4cout<<"*** Sponsored by: ESA/ESTEC contract contract 21435/08/NL/AT ***"<<std::endl;	
	G4cout<<"****************************************************************"<<std::endl;
	welcome_message=false;
  }	
  
  //Replace the user defined actions by the adjoint actions
  //---------------------------------------------------------
  SetAdjointPrimaryRunAndStackingActions();
  SetRestOfAdjointActions();
  
  //Update the list of primaries
  //-----------------------------
  theAdjointPrimaryGeneratorAction->UpdateListOfPrimaryParticles();

  adjoint_sim_mode=true;
  
  ID_of_last_particle_that_reach_the_ext_source=0;
  
  //Make the run
  //------------
  
  nb_evt_of_last_run =nb_evt;
  G4RunManager::GetRunManager()->BeamOn(theAdjointPrimaryGeneratorAction->GetNbOfAdjointPrimaryTypes()*2*nb_evt);

  //Restore the user defined actions
  //--------------------------------
  ResetRestOfUserActions(); 
  ResetUserPrimaryRunAndStackingActions();
  adjoint_sim_mode=false; 

  /*
  //Register the weight vector
  //--------------------------
  std::ofstream FileOutputElectronWeight("ElectronWeight.txt", std::ios::out);
  FileOutputElectronWeight<<std::setiosflags(std::ios::scientific);
  FileOutputElectronWeight<<std::setprecision(6);
  G4bool aBool = electron_last_weight_vector->Store(FileOutputElectronWeight, true);
  FileOutputElectronWeight.close();
  
  std::ofstream FileOutputProtonWeight("ProtonWeight.txt", std::ios::out);
  FileOutputProtonWeight<<std::setiosflags(std::ios::scientific);
  FileOutputProtonWeight<<std::setprecision(6);
  aBool = proton_last_weight_vector->Store(FileOutputProtonWeight, true);
  FileOutputProtonWeight.close();
  
  std::ofstream FileOutputGammaWeight("GammaWeight.txt", std::ios::out);
  FileOutputGammaWeight<<std::setiosflags(std::ios::scientific);
  FileOutputGammaWeight<<std::setprecision(6);
  aBool = gamma_last_weight_vector->Store(FileOutputGammaWeight, true);
  FileOutputGammaWeight.close();
  */  
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetRestOfAdjointActions()
{  
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  
  if (!user_action_already_defined) DefineUserActions();
 
 //Replace the user action by the adjoint actions
 //------------------------------------------------- 
  
  theRunManager->SetUserAction(theAdjointEventAction);
  theRunManager->SetUserAction(theAdjointSteppingAction);
  theRunManager->SetUserAction(theAdjointTrackingAction);  
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointPrimaryRunAndStackingActions()
{  
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  
  if (!user_action_already_defined) DefineUserActions();
  
 //Replace the user action by the adjoint actions
 //------------------------------------------------- 
  
  theRunManager->SetUserAction(theAdjointRunAction);
  theRunManager->SetUserAction(theAdjointPrimaryGeneratorAction);
  theRunManager->SetUserAction(theAdjointStackingAction);
  if (use_user_StackingAction)	theAdjointStackingAction->SetUserFwdStackingAction(fUserStackingAction);
  else theAdjointStackingAction->SetUserFwdStackingAction(0);
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::ResetRestOfUserActions()
{
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  
  //Restore the user defined actions
  //-------------------------------
 
  theRunManager->SetUserAction(fUserEventAction);
  theRunManager->SetUserAction(fUserSteppingAction);
  theRunManager->SetUserAction(fUserTrackingAction);
}

////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::ResetUserPrimaryRunAndStackingActions()
{ 
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  //Restore the user defined actions
  //-------------------------------
  theRunManager->SetUserAction(fUserRunAction);
  theRunManager->SetUserAction(fUserPrimaryGeneratorAction); 
  theRunManager->SetUserAction(fUserStackingAction);
}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::DefineUserActions()
{ 
   G4RunManager* theRunManager =  G4RunManager::GetRunManager();
   fUserTrackingAction= const_cast<G4UserTrackingAction* >( theRunManager->GetUserTrackingAction() );
   fUserEventAction= const_cast<G4UserEventAction* >( theRunManager->GetUserEventAction() ); 
   fUserSteppingAction= const_cast<G4UserSteppingAction* >( theRunManager->GetUserSteppingAction() );
   fUserPrimaryGeneratorAction= const_cast<G4VUserPrimaryGeneratorAction* >( theRunManager->GetUserPrimaryGeneratorAction() );
   fUserRunAction= const_cast<G4UserRunAction*>( theRunManager->GetUserRunAction() );
   fUserStackingAction= const_cast<G4UserStackingAction* >( theRunManager->GetUserStackingAction() );
   user_action_already_defined=true;	
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointTrackingMode(G4bool aBool)
{
  adjoint_tracking_mode = aBool;
 
  if (adjoint_tracking_mode) {
  	SetRestOfAdjointActions();
	theAdjointStackingAction->SetAdjointMode(true);
	theAdjointStackingAction->SetKillTracks(false);
	
  }	
  else {
  	
  	ResetRestOfUserActions();
	theAdjointStackingAction->SetAdjointMode(false);
	if (GetDidAdjParticleReachTheExtSource()){
		theAdjointStackingAction->SetKillTracks(false);
		RegisterAtEndOfAdjointTrack();
	}
	else theAdjointStackingAction->SetKillTracks(true);
  }	
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointSimManager::GetDidAdjParticleReachTheExtSource()
{
  return theAdjointSteppingAction->GetDidAdjParticleReachTheExtSource();
}
///////////////////////////////////////////////////////////////////////////////
//
std::vector<G4ParticleDefinition*>  G4AdjointSimManager::GetListOfPrimaryFwdParticles()
{
  return theAdjointPrimaryGeneratorAction->GetListOfPrimaryFwdParticles();
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::RegisterAtEndOfAdjointTrack()
{
  last_pos = theAdjointSteppingAction->GetLastPosition();  
  last_direction = theAdjointSteppingAction->GetLastMomentum();
  last_direction /=last_direction.mag();
  last_cos_th =  last_direction.z();
  G4ParticleDefinition* aPartDef= theAdjointSteppingAction->GetLastPartDef(); 
  	
  last_fwd_part_name= aPartDef->GetParticleName();
 
  last_fwd_part_name.remove(0,4);
  
  last_fwd_part_PDGEncoding=G4ParticleTable::GetParticleTable()->FindParticle(last_fwd_part_name)->GetPDGEncoding();
	
  std::vector<G4ParticleDefinition*> aList = theAdjointPrimaryGeneratorAction->GetListOfPrimaryFwdParticles();
  last_fwd_part_index=-1;
  size_t i=0;
  while(i<aList.size() && last_fwd_part_index<0) {
  	if (aList[i]->GetParticleName() == last_fwd_part_name) last_fwd_part_index=i;
	i++;
  }
  
  last_ekin = theAdjointSteppingAction->GetLastEkin();
  last_ekin_nuc = last_ekin;
  if (aPartDef->GetParticleType() == "adjoint_nucleus") {
  	nb_nuc=double(aPartDef->GetBaryonNumber());
  	last_ekin_nuc /=nb_nuc;
  }

  last_weight = theAdjointSteppingAction->GetLastWeight(); 
  
  /* G4PhysicsLogVector* theWeightVector=0;
  if (last_fwd_part_name =="e-")  theWeightVector=electron_last_weight_vector;
  else if (last_fwd_part_name =="gamma") theWeightVector=gamma_last_weight_vector;
  else if (last_fwd_part_name =="proton") theWeightVector=proton_last_weight_vector;
  
  if (theWeightVector){

	size_t ind =  size_t(std::log10(last_weight/theAdjointPrimaryWeight)*10. + 200);
   	G4double low_val =theWeightVector->GetLowEdgeEnergy(ind);
   	G4bool aBool = true;
   	G4double bin_weight = theWeightVector->GetValue(low_val, aBool)+1.;
   	theWeightVector->PutValue(ind, bin_weight);
  }
  */
  /*if ((last_weight/theAdjointPrimaryWeight)>1.) last_weight*=1000. ;
  else if ( (last_weight/theAdjointPrimaryWeight)>0.1) last_weight*=100. ;
  else if ( (last_weight/theAdjointPrimaryWeight)>0.01) last_weight*=10. ;*/
  
  
  //G4cout <<"Last Weight "<<last_weight<<'\t'<<theAdjointPrimaryWeight<<'\t'<<last_weight/theAdjointPrimaryWeight<<std::endl;
  /*if (last_weight/theAdjointPrimaryWeight >10.) {
  	G4cout<<"Warning a weight increase by a factor : "<<last_weight/theAdjointPrimaryWeight<<std::endl;
  }
  */
 
   ID_of_last_particle_that_reach_the_ext_source++;
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineSphericalExtSource(G4double radius, G4ThreeVector pos)
{	
   G4double area;
   return G4AdjointCrossSurfChecker::GetInstance()->AddaSphericalSurface("ExternalSource", radius, pos, area);
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineSphericalExtSourceWithCentreAtTheCentreOfAVolume(G4double radius, const G4String& volume_name)
{
   G4double area;
   G4ThreeVector center;
   return G4AdjointCrossSurfChecker::GetInstance()->AddaSphericalSurfaceWithCenterAtTheCenterOfAVolume( "ExternalSource", radius, volume_name,center, area);
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineExtSourceOnTheExtSurfaceOfAVolume(const G4String& volume_name)
{
   G4double area;
   return G4AdjointCrossSurfChecker::GetInstance()->AddanExtSurfaceOfAvolume( "ExternalSource", volume_name,area);
}
///////////////////////////////////////////////////////////////////////////////
//
void  G4AdjointSimManager::SetExtSourceEmax(G4double Emax)
{
  theAdjointSteppingAction->SetExtSourceEMax(Emax);
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineSphericalAdjointSource(G4double radius, G4ThreeVector pos)
{	
   G4double area;
   G4bool aBool = G4AdjointCrossSurfChecker::GetInstance()->AddaSphericalSurface("AdjointSource", radius, pos, area); 
   theAdjointPrimaryGeneratorAction->SetSphericalAdjointPrimarySource(radius, pos); 
   area_of_the_adjoint_source=area;
   return aBool; 	
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineSphericalAdjointSourceWithCentreAtTheCentreOfAVolume(G4double radius, const G4String& volume_name)
{
   G4double area;
   G4ThreeVector center;
   G4bool aBool = G4AdjointCrossSurfChecker::GetInstance()->AddaSphericalSurfaceWithCenterAtTheCenterOfAVolume( "AdjointSource", radius, volume_name,center, area);
   theAdjointPrimaryGeneratorAction->SetSphericalAdjointPrimarySource(radius, center);
   area_of_the_adjoint_source=area;
   return aBool;
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool  G4AdjointSimManager::DefineAdjointSourceOnTheExtSurfaceOfAVolume(const G4String& volume_name)
{
   G4double area;
   G4bool aBool = G4AdjointCrossSurfChecker::GetInstance()->AddanExtSurfaceOfAvolume( "AdjointSource", volume_name,area); 
   area_of_the_adjoint_source=area; 
   if (aBool) { 
   	theAdjointPrimaryGeneratorAction->SetAdjointPrimarySourceOnAnExtSurfaceOfAVolume(volume_name);
   }
   return aBool;
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointSourceEmin(G4double Emin)
{
  theAdjointPrimaryGeneratorAction->SetEmin(Emin);
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointSourceEmax(G4double Emax)
{
  theAdjointPrimaryGeneratorAction->SetEmax(Emax);
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::ConsiderParticleAsPrimary(const G4String& particle_name)
{
  theAdjointPrimaryGeneratorAction->ConsiderParticleAsPrimary(particle_name);
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::NeglectParticleAsPrimary(const G4String& particle_name)
{
  theAdjointPrimaryGeneratorAction->NeglectParticleAsPrimary(particle_name);
}
///////////////////////////////////////////////////////////////////////////////
//
/*void G4AdjointSimManager::SetPrimaryIon(G4int Z, G4int A)
{
  theAdjointPrimaryGeneratorAction->SetPrimaryIon(Z, A);
}
*/
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetPrimaryIon(G4ParticleDefinition* adjointIon, G4ParticleDefinition* fwdIon)
{
  theAdjointPrimaryGeneratorAction->SetPrimaryIon(adjointIon, fwdIon);
}
///////////////////////////////////////////////////////////////////////////////
//
const G4String& G4AdjointSimManager::GetPrimaryIonName()
{
  return theAdjointPrimaryGeneratorAction->GetPrimaryIonName();
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::RegisterAdjointPrimaryWeight(G4double aWeight)
{
  theAdjointPrimaryWeight = aWeight;
  theAdjointSteppingAction->SetPrimWeight(aWeight);
} 

///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointEventAction(G4UserEventAction* anAction)
{
  theAdjointEventAction = anAction;
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointSteppingAction(G4UserSteppingAction* anAction)
{
  theAdjointSteppingAction->SetUserAdjointSteppingAction(anAction);
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointStackingAction(G4UserStackingAction* anAction)
{
  theAdjointStackingAction->SetUserAdjointStackingAction(anAction);
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointTrackingAction(G4UserTrackingAction* anAction)
{
  theAdjointTrackingAction=anAction;
}
///////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSimManager::SetAdjointRunAction(G4UserRunAction* anAction)
{
  theAdjointRunAction=anAction;
} 
///////////////////////////////////////////////////////////////////////////////
//  
