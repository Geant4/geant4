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
// $Id: G4VNestedParameterisation.cc,v 1.3 2005-03-03 17:08:16 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VNestedParamterisation
//
// Class description:
//
// Parameterisation class, that can use parent volume information 
// eg to compute the transformation and the dimensions of parameterised 
// volumes, given a replication number and parent volume information.

// History:
// 04 Feb 05 J.Apostolakis Re-enabling the parameterisation using parent info
// --------------------------------------------------------------------
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VTouchable.hh"

#include "G4PhysicalTouchable.hh"

#include "G4VNestedParameterisation.hh" 

G4VPhysicalVolume* 
G4VNestedParameterisation::
ObtainParts( G4VPhysicalVolume* thisVolPlus,     // Input 
	     const G4String method, 
	     const G4VTouchable** pPtrTouchableParent)  const     // Returned
{
  G4PhysicalTouchable* pPhysTouchable=0;
  pPhysTouchable=  dynamic_cast<G4PhysicalTouchable *>(thisVolPlus); 
  
  if( ! pPhysTouchable ) 
    ReportErrorInTouchable( method, thisVolPlus ); 
  // ParameterisedNavigator must provide a G4PhysicalTouchable !!
  
  G4VPhysicalVolume* pCurrentVol= pPhysTouchable->GetCurrentVolume(); 
  // const G4VTouchable* 
  *pPtrTouchableParent= pPhysTouchable->GetParentTouchable();  // Parent Touch

  return pCurrentVol; 
}

const G4VPhysicalVolume* 
G4VNestedParameterisation::
ObtainParts( const G4VPhysicalVolume* thisVolPlus,     // Input 
	     const G4String method, 
	     const G4VTouchable** pPtrTouchableParent)  const     // Returned
{
  const G4PhysicalTouchable* pPhysTouchable=0;
  pPhysTouchable=  dynamic_cast<const G4PhysicalTouchable *>(thisVolPlus); 
  
  if( ! pPhysTouchable ) 
    ReportErrorInTouchable( method, thisVolPlus ); 
  // ParameterisedNavigator must provide a G4PhysicalTouchable !!
  
  const G4VPhysicalVolume* pCurrentVol= pPhysTouchable->GetCurrentVolume(); 
  // const G4VTouchable* 
  *pPtrTouchableParent= pPhysTouchable->GetParentTouchable();  // Parent Touch

  return pCurrentVol; 
}

G4Material* G4VNestedParameterisation::ComputeMaterial(const G4int no, 
						       G4VPhysicalVolume *currVolPlus)
{
  const G4String methodName("G4VNestedParameterisation::ComputeMaterial(int, physVol)");
  G4VPhysicalVolume* pCurrentVol= 0;  
  const G4VTouchable* pTouchableParent=0; 
  pCurrentVol= ObtainParts( currVolPlus, methodName, &pTouchableParent); 

  // Call the implementing virtual method
  return this->ComputeMaterial( no, pCurrentVol, pTouchableParent);  
}

void 
G4VNestedParameterisation::ReportErrorInTouchable(const G4String     method, 
						  const G4VPhysicalVolume* thisVol
) const
{
   const G4PhysicalTouchable* pPhysTouchable=0; 
   pPhysTouchable=  dynamic_cast<const G4PhysicalTouchable *>(thisVol); 

   G4cerr << "      Problem with cast of volume " << thisVol->GetName() 
	                  << " Address: " << thisVol << G4endl; 
   G4cerr << "       which results in phys-touch addr " << pPhysTouchable << G4endl;
   if( ! pPhysTouchable ){ 
     G4cerr << "ERROR - Illegal call to G4VNestedParameterisation method: ComputeDimensions, ComputeTransformation or ComputeMaterial()" << G4endl
	    << "        Physical volume is not 'dressed' to provide parent touchable. Cross cast failed!" << G4endl;
   
     G4Exception( method, "NotApplicable",
		  FatalException, "Failed cross cast: Illegal call, missing information.");
   }
}

G4VSolid*  G4VNestedParameterisation::ComputeSolid(const G4int, 
				     G4VPhysicalVolume  *pvol)    // currentVol
				     
{ 
  return pvol->GetLogicalVolume()->GetSolid(); 
}

