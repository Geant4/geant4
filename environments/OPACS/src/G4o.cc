//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4o.cc,v 1.7.4.1 2001/06/28 19:06:33 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
/*
#define DEBUG
*/

#include <stdio.h>
#include <string.h>

// Geant4 :
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4PhysicalVolumeModel.hh>
#include <G4ModelingParameters.hh>
#include <G4RunManager.hh>
#include <G4Trajectory.hh>

// For G4UI :
#include <G4UImanager.hh>
//#include <G4ParticleTblMessenger.hh>
#include <G4ParticleTable.hh>

// Co :
#include <CPrinter.h>
#include <CMemory.h>
#include <CString.h>
#include <CText.h>
#include <CFile.h>
#include <CList.h>
#include <OType.h>
#include <OShell.h>
#include <OProcess.h>

// Go :
#include <Go.h>
#include <GoTypes.h>

// G4o :
#include <G4oScene.hh>
#include <G4oDrawer.hh>
#include <G4o.h>

#include <G4oCommon.hh>

static void         ClearClass                    ();
static OIdentifier* GetPhysicalVolumeIdentifiers  (OType);
static int GetPhysicalVolumeProperty (OIdentifier,OType,OProperty,void*,int*);
static ONode        RepresentPhysicalVolume       (OIdentifier,OType);

static OIdentifier* GetTrajectoryIdentifiers      (OType);
static int GetTrajectoryProperty (OIdentifier,OType,OProperty,void*,int*);
static ONode        RepresentTrajectory           (OIdentifier,OType);


typedef unsigned long Ulong;

static struct {
  void**    extent;
  G4oScene* g4oScene;
  OShell    osh;
} Class = {NULL,NULL,NULL};
/***************************************************************************/
void ClearClass (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  CListDelete                          (Class.extent);
  Class.extent                         = NULL;
  if(Class.g4oScene!=NULL)             delete Class.g4oScene;
  Class.g4oScene                       = NULL;
  Class.osh                            = NULL;
}
/***************************************************************************/
void G4oAddCommands (
 void*    a_osh
)
/***************************************************************************/
/*
  Some G4 commands :
    G4> /particle/
    G4> /particle/dump
    G4> /particle/Verbose
  The definitions of these commands ar in G4o/src/G4oMessenger.cc code.
*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(OCommandGetIdentifier("G4")!=NULL) return; /*Done*/

  GoSetTypes ();
  SetTypes   ();

  if(OTypeGetIdentifier("PhysicalVolume")==NULL)
    {
      OType otype;
      otype                          = OTypeCreate ("PhysicalVolume");
      OTypeSetGetIdentifiersFunction (otype,GetPhysicalVolumeIdentifiers);
      OTypeSetGetPropertyFunction    (otype,GetPhysicalVolumeProperty);
      OTypeSetClearClassFunction     (otype,ClearClass);
      OTypeSetRepresentFunction      (otype,(OTypeRepresentFunction)RepresentPhysicalVolume);
      OTypeAddNewProperty            (otype,"identifier" ,OPropertyUnsignedLong    ,NULL);
      OTypeAddNewProperty            (otype,"name"       ,OPropertyString          ,NULL);
      
      otype                          = OTypeCreate ("Trajectory");
      OTypeSetGetIdentifiersFunction (otype,GetTrajectoryIdentifiers);
      OTypeSetGetPropertyFunction    (otype,GetTrajectoryProperty);
      OTypeSetRepresentFunction      (otype,
				(OTypeRepresentFunction)RepresentTrajectory);
      OTypeAddNewProperty (otype,"identifier"  ,OPropertyUnsignedLong,NULL);
      OTypeAddNewProperty (otype,"trackID"     ,OPropertyInteger     ,NULL);
      OTypeAddNewProperty (otype,"parentID"    ,OPropertyInteger     ,NULL);
      OTypeAddNewProperty (otype,"particleName",OPropertyString      ,NULL);
      OTypeAddNewProperty (otype,"charge"      ,OPropertyDouble      ,NULL);
    }


  OShellAddNewCommand ((OShell)a_osh,"G4o/G4",Execute_G4);

  Class.osh = (OShell)a_osh;  
}
/***************************************************************************/
void G4oExecuteScript (
 char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  ExecuteScript (a_string);
}
/***************************************************************************/
OShell G4oGetShell (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return Class.osh;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
OIdentifier* GetPhysicalVolumeIdentifiers (
 OType a_type
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  int objn = G4PhysicalVolumeStore::GetInstance()->entries();
#ifdef DEBUG
  printf ("debug: GetPhysicalVolumeIdentifiers : %d objects.\n",objn);
#endif
  if(objn<=0)        return NULL;
  CListDelete        (Class.extent);
  Class.extent       = CListCreate(objn);
  if(Class.extent==NULL)  return NULL;
  for(int count=0;count<objn;count++) {
    Class.extent[count] = (OIdentifier)(count+1);
  }
  a_type             = NULL;
  return             Class.extent;
}
/***************************************************************************/
int GetPhysicalVolumeProperty (
 OIdentifier a_obj
,OType  This
,OProperty  a_prop
,void*  a_addr 
,int*   a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_number!=NULL) *a_number = 0;
  if(This==NULL) return 0;
  char* name = OPropertyGetName (a_prop);
  if(name==NULL) return 0;

  if(strcmp(name,"identifier")==0) {
    *((Ulong*)a_addr) = (Ulong)a_obj;
  } else if(strcmp(name,"name")==0) {
    *((char**)a_addr) = CStringDuplicate((char*)(G4PhysicalVolumeStore::GetInstance()->at((size_t)a_obj-1)->GetName().data()));
    return            FREE_BLOCK;
  } else {
    return 0;
  }
  return 1;
}
/***************************************************************************/
ONode RepresentPhysicalVolume (
 OIdentifier a_obj
,OType This
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(This==NULL)     return NULL;

  G4VPhysicalVolume* physicalVolume = 
    G4PhysicalVolumeStore::GetInstance()->at((size_t)a_obj-1);

  if(Class.g4oScene==NULL) {
    Class.g4oScene = new G4oScene;
    if(Class.g4oScene==NULL) return NULL;
  }

  char* string = CStringCreateF (15+64,"PhysicalVolume/%lu",a_obj);
  Class.g4oScene->SetNodeName (string);
  CStringDelete               (string);

  G4ModelingParameters mParams (0,
				G4ModelingParameters::wireframe,
				true,  // Global culling.
				true,  // Cull invisible volumes.
				false, // Density culling.
				0.,    // Density (not relevant if density 
				       // culling false).
				true,  // Cull daughters of opaque mothers.
				24);   // No of sides (not relevant for 
                                       // this operation).

  G4PhysicalVolumeModel model  (physicalVolume,
				G4PhysicalVolumeModel::UNLIMITED,
				G4Transform3D::Identity,
				&mParams);
  model.DescribeYourselfTo     (*Class.g4oScene);

  return                       Class.g4oScene->GetNode();
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
OIdentifier* GetTrajectoryIdentifiers (
 OType a_type
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  const G4Event* event  = G4RunManager::GetRunManager() -> GetCurrentEvent();
  if(event==NULL) return 0;
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  if(trajectoryContainer==NULL) return 0;
  int objn = trajectoryContainer->entries();
#ifdef DEBUG
  printf ("debug: GetTrajectoryIdentifiers : event : %lu, %d objects.\n",
	  event,objn);
#endif
  if(objn<=0)        return NULL;
  CListDelete        (Class.extent);
  Class.extent       = CListCreate(objn);
  if(Class.extent==NULL)  return NULL;
  for(int count=0;count<objn;count++) {
    Class.extent[count] = (OIdentifier)(count+1);
  }
  a_type             = NULL;
  return             Class.extent;
}
/***************************************************************************/
int GetTrajectoryProperty (
 OIdentifier a_obj
,OType  This
,OProperty  a_prop
,void*  a_addr 
,int*   a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_number!=NULL)            *a_number = 0;
  if(This==NULL)                return 0;
  char* name = OPropertyGetName (a_prop);
  if(name==NULL)                return 0;
  const G4Event* event  = G4RunManager::GetRunManager() -> GetCurrentEvent();
  if(event==NULL)               return 0;
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  if(trajectoryContainer==NULL) return 0;

  G4Trajectory* trajectory = 
    (G4Trajectory*)((*trajectoryContainer)[(size_t)a_obj-1]);

  if(strcmp(name,"identifier")==0) {
    *((Ulong*)a_addr) = (Ulong)a_obj;
  } else if(strcmp(name,"trackID")==0) {
    *((int*)a_addr)   = trajectory->GetTrackID();
  } else if(strcmp(name,"parentID")==0) {
    *((int*)a_addr)   = trajectory->GetParentID();
  } else if(strcmp(name,"particleName")==0) {
    *((char**)a_addr) = CStringDuplicate(
             (char*)(trajectory->GetParticleName().data()));
    return            FREE_BLOCK;
  } else if(strcmp(name,"charge")==0) {
    *((double*)a_addr) = trajectory->GetCharge();
  } else {
    return 0;
  }
  return 1;
}
/***************************************************************************/
ONode RepresentTrajectory (
 OIdentifier a_obj
,OType This
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(This==NULL)                return NULL;
  const G4Event* event  = G4RunManager::GetRunManager() -> GetCurrentEvent();
  if(event==NULL)               return NULL;
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  if(trajectoryContainer==NULL) return NULL;
  G4oDrawer*  drawer = G4oDrawer::GetInstance();
  if(drawer==NULL) return NULL;

#ifdef DEBUG
  printf ("debug: RepresentTrajectory : %d.\n",a_obj);
#endif

  ONode node = ONodeCreateF (11+64,"Trajectory/%lu",a_obj);

  drawer->SetNode (node);
  G4Trajectory* trajectory = 
    (G4Trajectory*)((*trajectoryContainer)[(size_t)a_obj-1]);

  trajectory->DrawTrajectory(0); 
  
  return          node;
}
#include <G4oCommon.icc>







