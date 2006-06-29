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
// Satoshi Tanaka  31th May 2001
// A scene handler to dump geometry hierarchy to GAG.

#include "G4Types.hh"
#include <sstream>

#include "G4GAGTreeSceneHandler.hh"

#include "G4GAGTree.hh"
#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"

// #define DEBUG_GAG_TREE

const G4String  GAG_TREE_WORLD_VOLUME_NAME ( "WORLD" ); 
const G4String  GAG_TREE_BEGIN_DTREE       ( "@@DTREEBegin" ); 
const G4String  GAG_TREE_END_DTREE         ( "@@DTREEEnd" ); 


//-----
G4GAGTreeSceneHandler::G4GAGTreeSceneHandler(G4VGraphicsSystem& system,
						 const G4String& name):
  G4VTreeSceneHandler(system, name) 
{
  ClearPVList();
}

//-----
G4GAGTreeSceneHandler::~G4GAGTreeSceneHandler () 
{
  ClearPVList();
}

//-----
void G4GAGTreeSceneHandler::ClearPVList(void)
{ 
  fPVNameList.Clear() ;  

  fPrevDepth = 0 ; 
  fPrevAbsPVName = "";
  fPVCounter = 0 ;
} 

//-----
void G4GAGTreeSceneHandler::InitializePVList(void)
{ 
  fPVNameList.Clear(); 
  fPVNameList.Push( GAG_TREE_WORLD_VOLUME_NAME ); 

  fPrevDepth = 0 ; 
  fPrevAbsPVName = ""; 
  fPVCounter = 0 ;
} 

//-----
void G4GAGTreeSceneHandler::BeginModeling () 
{
  G4VTreeSceneHandler::BeginModeling ();  

  // Initialize the stack the the tree generation
  InitializePVList();

  // Get the current verbosity level 
  const G4GAGTree* pSystem = (G4GAGTree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();

  // Beginning of the tree output
  G4cout << "\nG4GAGTreeSceneHandler::BeginModeling:"
    "\n  set verbosity with \"/vis/GAGTree/verbose <verbosity>\":"
    "\n  <  10: - does not print daughters of repeated placements."
    "\n         - does not repeat replicas."
    "\n  >= 10: prints all physical volumes."
    "\n          FORMAT: PV_name.copy_no.index"
    "\n  For level of detail add:"
    "\n  >=  0: prints physical volume name and copy number."
    "\n          FORMAT: PV_name.copy_no.index"
    "\n  >=  1: prints logical volume name."
    "\n          FORMAT: PV_name.copy_no.LV_name.index"
    "\n  >=  2: prints solid name and type."
    "\n          FORMAT: PV_name.copy_no.LV_name.solid_name.solid_type.index"
    "\n  Note: all culling, if any, is switched off so all volumes are seen."
    "\n  Now printing with verbosity " << verbosity << G4endl;

  // Declare to start GAG tree
  G4cout << GAG_TREE_BEGIN_DTREE << G4endl;

}

void G4GAGTreeSceneHandler::EndModeling () 
{
  // Dump the last PV node
  if( fPrevAbsPVName != "" ) { G4cout << fPrevAbsPVName << G4endl; }

  // Declare to end GAG tree
  G4cout << GAG_TREE_END_DTREE << G4endl;
  G4cout << "G4GAGTreeSceneHandler::EndModeling" << G4endl;

  // Clear the sets
  fLVSet.clear();
  fReplicaSet.clear();

  // Clear the stack for the tree generation
  ClearPVList();

  G4VTreeSceneHandler::EndModeling ();  
}

void G4GAGTreeSceneHandler::RequestPrimitives (const G4VSolid&) 
{
  // Protection if solid drawn not from physical volume tree...(JA)
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) return;
  G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  G4int currentDepth = pPVModel->GetCurrentDepth();

//////////////////////////////
  G4String        cur_abs_pv_name ;
  G4String        pv_name_tmp ;
  G4String        cur_pv_name ( pCurrentPV->GetName() ) ; 
//////////////////////////////

  const G4GAGTree* pSystem = (G4GAGTree*)GetGraphicsSystem();
  const G4int verbosity = pSystem->GetVerbosity();
  const G4int detail = verbosity % 10;

  if (verbosity < 10 && fReplicaSet.find(pCurrentPV) != fReplicaSet.end()) {
    // Ignore if an already treated replica.
    pPVModel->CurtailDescent();
    return;
  }


///////////////////////////////////////////////////////////
  // Add the extension ".details.index" to the current PV node

    // Step 1: Initialize the extension
  std::ostringstream ost;

    // Step 2: Generate the extension
      // copy number 
  ost << "." << pCurrentPV->GetCopyNo() ;
  if (detail >= 1) {
      // LV name
    ost << "." << pCurrentLV->GetName() ; 
  }
  if (detail >= 2) {
      // Solid info
    ost << "." << pCurrentLV->GetSolid()->GetName();
    ost	<< "." << pCurrentLV->GetSolid()->GetEntityType() ;
  }
      // Tree index (0 for the world) 
  ost << "." << fPVCounter; 

    // Step 3: Add the extention to the current PV name
  cur_pv_name += ost.str() ;            

  // End of adding extension


  // Search a mother PV node for the current PV
  //  in the direction to the root node.
  //  depth_mother = depth_current - 1 
  if( currentDepth >  fPrevDepth )
  {
    // Do nothing (The mother is the previous PV node)

  } else { 

    while (1) 
    {
      // Delete the head item of the list
      fPVNameList.Pop();

      // Is the right mother at the head of the list?
      //  Note: The mother of the world has depth -1.
      //        It happens when the current node is the world volume, 
      //        and so the list becomes empty with the popping above.

      G4int trial_mother_depth = fPVNameList.GetHeadIndex() ;
      if ( currentDepth > trial_mother_depth ) break ;
    } 

  } // if-else

  // Add the curent PV to the tree as a new node/leaf
  fPVNameList.Push( cur_pv_name );

  // Make the absolute path name
  fPVNameList.ToTail() ;
  while ( fPVNameList.GetItem ( pv_name_tmp ) ) {
        cur_abs_pv_name +=  "/" ;
	cur_abs_pv_name += pv_name_tmp ;
        fPVNameList.Upward();
  }


  // Add "/"  to the PREVIOUS PV node name if it is not a leaf.
  if( fPrevAbsPVName != "" ) 
  { 
    // Add "/"
    if( currentDepth > fPrevDepth ) {fPrevAbsPVName += "/" ;}
  } 

  // Print the PREVIOUS node after node/leaf distiction 
  //  is made clear by reading the current node.
   if( fPrevAbsPVName != "" ) { G4cout << fPrevAbsPVName ; }
///////////////////////////////////////////////////////////


  if (pCurrentPV->IsReplicated()) {
    fReplicaSet.insert(pCurrentPV);  // Record new replica volume.
    if (verbosity < 10) {
      // Add printing for replicas (when replicas are ignored)...
      EAxis axis;
      G4int nReplicas;
      G4double width;
      G4double offset;
      G4bool consuming;
      pCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
      G4VPVParameterisation* pP = pCurrentPV->GetParameterisation();
#if defined DEBUG_GAG_TREE
      G4cout << "_(" << nReplicas;
#endif
      if (pP) {
#if defined DEBUG_GAG_TREE
	G4cout << " parametrised volumes)";
#endif
      }
      else {
#if defined DEBUG_GAG_TREE
	G4cout << " replicas)";
#endif
      }
    }
  }
  else {
    if (fLVSet.find(pCurrentLV) != fLVSet.end()) {
      if (verbosity <  10) {
	// Add printing for repeated placement...
#if defined DEBUG_GAG_TREE
	G4cout << " (repeated placement)";
#endif
	// Ignore if an already treated logical volume.
	pPVModel->CurtailDescent();
	G4cout << G4endl;
	return;
      }
    }
  }

  if (fLVSet.find(pCurrentLV) == fLVSet.end()) {
    fLVSet.insert(pCurrentLV);  // Record new logical volume.
  }


////////////////////////
  // End of printing
  if( fPrevAbsPVName != "" ) { G4cout << G4endl;}

  // Prepare for the next call, i.e. the next PV.
  fPrevDepth     = currentDepth ;
  fPrevAbsPVName = cur_abs_pv_name ;
  fPVCounter++ ;
////////////////////////


  return;
}
