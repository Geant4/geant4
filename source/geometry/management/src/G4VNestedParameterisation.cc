//
// $Id: G4VNestedParameterisation.cc,v 1.1 2005-02-15 17:17:26 japost Exp $
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
#include "G4VTouchable.hh"

#include "G4PhysicalTouchable.hh"

#include "G4VNestedParameterisation.hh" 


// Implement standard PVParameterisation methods, 
//   in terms of the new methods with 'parent' information

G4VSolid*  
G4VNestedParameterisation::ComputeSolid(const G4int no, G4VPhysicalVolume *thisVol)
{
  G4PhysicalTouchable* pPhysTouchable;
  pPhysTouchable=  cross_cast<G4PhysicalTouchable*> thisVol; 
  
  if( ! pPhysTouchable ) 
    G4Exception(" ") ; // FATAL error -- 
  // ParameterisedNavigator must provide a G4PhysicalTouchable
  
  G4VPhysicalVolume pCurrentVol; 
  const G4VTouchable* pTouchableParent; 
  // currentVol= 
  this->ComputeSolid( no, pCurrentVol, pTouchableParent); 
}


void G4VNestedParameterisation::ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume *currPhysTouch )
{
  const G4String method("G4VNestedParameterisation::Transformation(int, pv)");

  G4PhysicalTouchable* pPhysTouchable;
  pPhysTouchable=  cross_cast<G4PhysicalTouchable*> thisVol; 
  
  if( ! pPhysTouchable ) 
    ReportErrorInTouchable(method); 
  // ParameterisedNavigator must provide a G4PhysicalTouchable
  
  G4VPhysicalVolume pCurrentVol; 
  const G4VTouchable* pTouchableParent; 
  // currentVol= 
  this->ComputeTransformation( no, pCurrentVol, pTouchableParent);   
}

G4Material* G4VNestedParameterisation::ComputeMaterial(const G4int, G4VPhysicalVolume *)
{
  const G4String Method("G4VNestedParameterisation::ComputeDimensions()");
}

void G4VNestedParameterisation::ReportErrorInTouchable(G4String method)
{
   G4cerr << "ERROR - Illegal call to G4VNestedParameterisation method: ComputeDimensions, ComputeTransformation or ComputeMaterial()" << G4endl
           << "        Physical volume is not 'dressed' to provide parent touchable. Cross cast failed!" << G4endl;
   G4Exception( method, "NotApplicable",
                FatalException, "Failed cross cast: Illegal call, missing information.");
}
