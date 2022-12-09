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
// class G4RegularNavigation implementation
//
// Author: Pedro Arce, May 2007
//
// --------------------------------------------------------------------

#include "G4RegularNavigation.hh"
#include "G4TouchableHistory.hh"
#include "G4PhantomParameterisation.hh"
#include "G4Material.hh"
#include "G4NormalNavigation.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"
#include "G4RegularNavigationHelper.hh"

//------------------------------------------------------------------
G4RegularNavigation::G4RegularNavigation()
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fMinStep = 101*kCarTolerance;
}


//------------------------------------------------------------------
G4RegularNavigation::~G4RegularNavigation()
{
}


//------------------------------------------------------------------
G4double G4RegularNavigation::
                    ComputeStep(const G4ThreeVector& localPoint,
                                const G4ThreeVector& localDirection,
                                const G4double currentProposedStepLength,
                                      G4double& newSafety,
                                      G4NavigationHistory& history,
                                      G4bool& validExitNormal,
                                      G4ThreeVector& exitNormal,
                                      G4bool& exiting,
                                      G4bool& entering,
                                      G4VPhysicalVolume *(*pBlockedPhysical),
                                      G4int& blockedReplicaNo)
{
  // This method is never called because to be called the daughter has to be
  // a regular structure. This would only happen if the track is in the mother
  // of voxels volume. But the voxels fill completely their mother, so when a
  // track enters the mother it automatically enters a voxel. Only precision
  // problems would make this method to be called

  G4ThreeVector globalPoint =
    history.GetTopTransform().InverseTransformPoint(localPoint);
  G4ThreeVector globalDirection =
         history.GetTopTransform().InverseTransformAxis(localDirection);

  G4ThreeVector localPoint2 = localPoint; // take away constantness

  LevelLocate( history, *pBlockedPhysical, blockedReplicaNo, 
               globalPoint, &globalDirection, true, localPoint2 );


  // Get in which voxel it is
  //
  G4VPhysicalVolume *motherPhysical, *daughterPhysical;
  G4LogicalVolume *motherLogical;
  motherPhysical = history.GetTopVolume();
  motherLogical = motherPhysical->GetLogicalVolume();
  daughterPhysical = motherLogical->GetDaughter(0);

  G4PhantomParameterisation * daughterParam =
        (G4PhantomParameterisation*)(daughterPhysical->GetParameterisation());
  G4int copyNo = daughterParam ->GetReplicaNo(localPoint,localDirection);

  G4ThreeVector voxelTranslation = daughterParam->GetTranslation( copyNo );
  G4ThreeVector daughterPoint = localPoint - voxelTranslation;


  // Compute step in voxel
  //
  return fnormalNav->ComputeStep(daughterPoint,
                                 localDirection,
                                 currentProposedStepLength,
                                 newSafety,
                                 history,
                                 validExitNormal,
                                 exitNormal,
                                 exiting,
                                 entering,
                                 pBlockedPhysical,
                                 blockedReplicaNo);
}


//------------------------------------------------------------------
G4double G4RegularNavigation::ComputeStepSkippingEqualMaterials(
                                      G4ThreeVector& localPoint,
                                const G4ThreeVector& localDirection,
                                const G4double currentProposedStepLength,
                                G4double& newSafety,
                                G4NavigationHistory& history,
                                G4bool& validExitNormal,
                                G4ThreeVector& exitNormal,
                                G4bool& exiting,
                                G4bool& entering,
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int& blockedReplicaNo,
                                G4VPhysicalVolume* pCurrentPhysical)
{
  G4RegularNavigationHelper::Instance()->ClearStepLengths();

  G4PhantomParameterisation *param =
    (G4PhantomParameterisation*)(pCurrentPhysical->GetParameterisation());

  if( !param->SkipEqualMaterials() )
  {
    return fnormalNav->ComputeStep(localPoint,
                                   localDirection,
                                   currentProposedStepLength,
                                   newSafety,
                                   history,
                                   validExitNormal,
                                   exitNormal,
                                   exiting,
                                   entering,
                                   pBlockedPhysical,
                                   blockedReplicaNo);
  }


  G4double ourStep = 0.;

  // To get replica No: transform local point to the reference system of the
  // param container volume
  //
  G4int ide = (G4int)history.GetDepth();
  G4ThreeVector containerPoint = history.GetTransform(ide)
                                 .InverseTransformPoint(localPoint);

  // Point in global frame
  //
  containerPoint = history.GetTransform(ide).InverseTransformPoint(localPoint);

  // Point in voxel parent volume frame
  //
  containerPoint = history.GetTransform(ide-1).TransformPoint(containerPoint);

  // Store previous voxel translation to move localPoint by the difference
  // with the new one
  //
  G4ThreeVector prevVoxelTranslation = containerPoint - localPoint;

  // Do not use the expression below: There are cases where the
  // fLastLocatedPointLocal does not give the correct answer
  // (particle reaching a wall and bounced back, particle travelling through
  // the wall that is deviated in an step, ...; these are pathological cases
  // that give wrong answers in G4PhantomParameterisation::GetReplicaNo()
  //
  // G4ThreeVector prevVoxelTranslation = param->GetTranslation( copyNo );

  G4int copyNo = param->GetReplicaNo(containerPoint,localDirection);

  G4Material* currentMate = param->ComputeMaterial( copyNo, nullptr, nullptr );
  G4VSolid* voxelBox = pCurrentPhysical->GetLogicalVolume()->GetSolid();

  G4VSolid* containerSolid = param->GetContainerSolid();
  G4Material* nextMate;
  G4bool bFirstStep = true;
  G4double newStep;
  G4double totalNewStep = 0.;

  // Loop while same material is found 
  //
  //
  fNumberZeroSteps = 0;
  for( G4int ii = 0; ii < fNoStepsAllowed+1; ++ii )
  {
    if( ii == fNoStepsAllowed ) {
      // Must kill this stuck track
      //
      G4ThreeVector pGlobalpoint = history.GetTransform(ide)
                                   .InverseTransformPoint(localPoint);
      std::ostringstream message;
      message << "G4RegularNavigation::ComputeStepSkippingEqualMaterials()"
        << "Stuck Track: potential geometry or navigation problem."
        << G4endl
        << "        Track stuck, moving for more than " 
        << ii << " steps" << G4endl
        << "- at point " << pGlobalpoint << G4endl
        << "        local direction: " << localDirection << G4endl;
      G4Exception("G4RegularNavigation::ComputeStepSkippingEqualMaterials()",
      "GeomRegNav1001",
      EventMustBeAborted,
      message);
    }
    newStep = voxelBox->DistanceToOut( localPoint, localDirection );
    fLastStepWasZero = (newStep<fMinStep);
    if( fLastStepWasZero )
    {
      ++fNumberZeroSteps;
#ifdef G4DEBUG_NAVIGATION
      if( fNumberZeroSteps > 1 )
      {
        G4ThreeVector pGlobalpoint = history.GetTransform(ide)
                                     .InverseTransformPoint(localPoint);
      std::ostringstream message;
      message.precision(16);
      message << "G4RegularNavigation::ComputeStepSkippingEqualMaterials(): another 'zero' step, # "
            << fNumberZeroSteps
            << ", at " << pGlobalpoint
            << ", nav-comp-step calls # " << ii
            << ", Step= " << newStep;
            G4Exception("G4RegularNavigation::ComputeStepSkippingEqualMaterials()",
                        "GeomRegNav1002", JustWarning, message,
                        "Potential overlap in geometry!");
      }
#endif
      if( fNumberZeroSteps > fActionThreshold_NoZeroSteps-1 )
      {
        // Act to recover this stuck track. Pushing it along direction
        //
        newStep = std::min(101*kCarTolerance*std::pow(10,fNumberZeroSteps-2),0.1);
#ifdef G4DEBUG_NAVIGATION
        G4ThreeVector pGlobalpoint = history.GetTransform(ide)
                                       .InverseTransformPoint(localPoint);
        std::ostringstream message;
        message.precision(16);
        message << "Track stuck or not moving." << G4endl
                      << "          Track stuck, not moving for " 
                      << fNumberZeroSteps << " steps" << G4endl
                      << "- at point " << pGlobalpoint
                      << " (local point " << localPoint << ")" << G4endl
                      << "        local direction: " << localDirection 
                      << "          Potential geometry or navigation problem !"
                      << G4endl
                      << "          Trying pushing it of " << newStep << " mm ...";
              G4Exception("G4RegularNavigation::ComputeStepSkippingEqualMaterials()",
                          "GeomRegNav1003", JustWarning, message,
                          "Potential overlap in geometry!");
#endif
      }
      if( fNumberZeroSteps > fAbandonThreshold_NoZeroSteps-1 )
      {
        // Must kill this stuck track
        //
        G4ThreeVector pGlobalpoint = history.GetTransform(ide)
                                          .InverseTransformPoint(localPoint);
        std::ostringstream message;
        message << "G4RegularNavigation::ComputeStepSkippingEqualMaterials()"
          << "Stuck Track: potential geometry or navigation problem."
          << G4endl
          << "        Track stuck, not moving for " 
          << fNumberZeroSteps << " steps" << G4endl
          << "- at point " << pGlobalpoint << G4endl	
          << "        local direction: " << localDirection << G4endl;
        G4Exception("G4RegularNavigation::ComputeStepSkippingEqualMaterials()",
              "GeomRegNav1004",
              EventMustBeAborted,
              message);
      }
    }
    else
    {
      // reset the zero step counter when a non-zero step was performed
      fNumberZeroSteps = 0;
    }
    if( (bFirstStep) && (newStep < currentProposedStepLength) )
    {
      exiting  = true;
    }
    bFirstStep = false;
 
    newStep += kCarTolerance;   // Avoid precision problems
    ourStep += newStep;
    totalNewStep += newStep;

    // Physical process is limiting the step, don't continue
    //
    if(std::fabs(totalNewStep-currentProposedStepLength) < kCarTolerance)
    { 
      return currentProposedStepLength;
    }
    if(totalNewStep > currentProposedStepLength) 
    { 
      G4RegularNavigationHelper::Instance()->
        AddStepLength(copyNo, newStep-totalNewStep+currentProposedStepLength);
      return currentProposedStepLength;
    }
    else
    {
      G4RegularNavigationHelper::Instance()->AddStepLength( copyNo, newStep );
    }


    // Move container point until wall of voxel
    //
    containerPoint += newStep*localDirection;
    if( containerSolid->Inside( containerPoint ) != kInside )
    {
      break;
    }

    // Get copyNo and translation of new voxel
    //
    copyNo = param->GetReplicaNo(containerPoint, localDirection);
    G4ThreeVector voxelTranslation = param->GetTranslation( copyNo );

    // Move local point until wall of voxel and then put it in the new voxel
    // local coordinates
    //
    localPoint += newStep*localDirection;
    localPoint += prevVoxelTranslation - voxelTranslation;

    prevVoxelTranslation = voxelTranslation;

    // Check if material of next voxel is the same as that of the current voxel
    nextMate = param->ComputeMaterial( copyNo, nullptr, nullptr );

    if( currentMate != nextMate ) { break; }
  }

  return ourStep;
}


//------------------------------------------------------------------
G4double
G4RegularNavigation::ComputeSafety(const G4ThreeVector& localPoint,
                                   const G4NavigationHistory& history,
                                   const G4double pMaxLength)
{
  // This method is never called because to be called the daughter has to be a
  // regular structure. This would only happen if the track is in the mother of
  // voxels volume. But the voxels fill completely their mother, so when a
  // track enters the mother it automatically enters a voxel. Only precision
  // problems would make this method to be called

  // Compute step in voxel
  //
  return fnormalNav->ComputeSafety(localPoint,
                                   history,
                                   pMaxLength );
}


//------------------------------------------------------------------
G4bool
G4RegularNavigation::LevelLocate( G4NavigationHistory& history,
                                  const G4VPhysicalVolume* ,
                                  const G4int ,
                                  const G4ThreeVector& globalPoint,
                                  const G4ThreeVector* globalDirection,
                                  const G4bool, // pLocatedOnEdge, 
                                  G4ThreeVector& localPoint )
{
  G4VPhysicalVolume *motherPhysical, *pPhysical;
  G4PhantomParameterisation *pParam;
  G4LogicalVolume *motherLogical;
  G4ThreeVector localDir;
  G4int replicaNo;
  
  motherPhysical = history.GetTopVolume();
  motherLogical = motherPhysical->GetLogicalVolume();
  
  pPhysical = motherLogical->GetDaughter(0);
  pParam = (G4PhantomParameterisation*)(pPhysical->GetParameterisation());
  
  // Save parent history in touchable history
  // ... for use as parent t-h in ComputeMaterial method of param
  //
  G4TouchableHistory parentTouchable( history ); 
  
  // Get local direction
  //
  if( globalDirection )
  {
    localDir = history.GetTopTransform().TransformAxis(*globalDirection);
  }
  else
  {
    localDir = G4ThreeVector(0.,0.,0.);
  }

  // Enter this daughter
  //
  replicaNo = pParam->GetReplicaNo( localPoint, localDir );

  if( replicaNo < 0 || replicaNo >= G4int(pParam->GetNoVoxels()) )
  {
    return false;
  }

  // Set the correct copy number in physical
  //
  pPhysical->SetCopyNo(replicaNo);
  pParam->ComputeTransformation(replicaNo,pPhysical);

  history.NewLevel(pPhysical, kParameterised, replicaNo );
  localPoint = history.GetTopTransform().TransformPoint(globalPoint);

  // Set the correct solid and material in Logical Volume
  //
  G4LogicalVolume *pLogical = pPhysical->GetLogicalVolume();
      
  pLogical->UpdateMaterial(pParam->ComputeMaterial(replicaNo,
                           pPhysical, &parentTouchable) );
  return true;
}
