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
//
// $Id$
//
//
// --------------------------------------------------------------------
//      GEANT 4 class implementation file
// --------------------------------------------------------------------

#include "G4ErrorPropagationNavigator.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorSurfaceTarget.hh"

//-------------------------------------------------------------------

G4ErrorPropagationNavigator::G4ErrorPropagationNavigator()
  : G4Navigator()
{
}

//-------------------------------------------------------------------

G4ErrorPropagationNavigator::~G4ErrorPropagationNavigator()
{
}

//-------------------------------------------------------------------

G4double G4ErrorPropagationNavigator::
ComputeStep ( const G4ThreeVector &pGlobalPoint,
              const G4ThreeVector &pDirection,
              const G4double pCurrentProposedStepLength,
                    G4double &pNewSafety )
{
  G4double Step = G4Navigator::ComputeStep(pGlobalPoint, pDirection,
                                           pCurrentProposedStepLength,
                                           pNewSafety);
  
  G4ErrorPropagatorData * g4edata
    = G4ErrorPropagatorData::GetErrorPropagatorData();

  if (g4edata !=0)
  {
    const G4ErrorTarget* target = g4edata->GetTarget();
    if( target != 0 )
    {
      G4double StepPlane= target->GetDistanceFromPoint(pGlobalPoint,pDirection);

      if( StepPlane < 0. ) // Negative means target is crossed, will not be found
      {
        StepPlane = DBL_MAX; 
      }
#ifdef G4VERBOSE
      if( G4ErrorPropagatorData::verbose() >= 4 )
      {
        G4cout << "G4ErrorPropagationNavigator::ComputeStep()" << G4endl
               << "  Target step: " << StepPlane
               << ", Transportation step: " << Step << G4endl;
        target->Dump( "G4ErrorPropagationNavigator::ComputeStep Target " );
      }
#endif

      if(StepPlane<Step)
      {
#ifdef G4VERBOSE
        if( G4ErrorPropagatorData::verbose() >= 2 )
        {
          G4cout << "G4ErrorPropagationNavigator::ComputeStep()" << G4endl
                 << "  TargetCloserThanBoundary: " << StepPlane << " < "
                 << Step << G4endl;
        }
#endif
        Step = StepPlane;
        g4edata->SetState(G4ErrorState_TargetCloserThanBoundary);
      }
      else
      {
        g4edata->SetState(G4ErrorState_Propagating);
      }
    }
  }
  pNewSafety = ComputeSafety(pGlobalPoint, pCurrentProposedStepLength);

#ifdef G4VERBOSE
  if( G4ErrorPropagatorData::verbose() >= 3 )
  {
    G4cout << "G4ErrorPropagationNavigator::ComputeStep()" << G4endl
           << "  Step: " << Step << ", ComputeSafety: " << pNewSafety
           << G4endl;
  }
#endif

  return Step;
}

//-------------------------------------------------------------------

G4double G4ErrorPropagationNavigator::
ComputeSafety( const G4ThreeVector &pGlobalpoint,
               const G4double pMaxLength,
               const G4bool keepState )
{
  G4double newSafety = G4Navigator::ComputeSafety(pGlobalpoint,
                                                  pMaxLength, keepState);

  G4ErrorPropagatorData *g4edata
    = G4ErrorPropagatorData::GetErrorPropagatorData();

  if (g4edata !=0)
  {
    const G4ErrorTarget* target = g4edata->GetTarget();
    if( target != 0 )
    {
      G4double distance = target->GetDistanceFromPoint(pGlobalpoint);
      
      if(distance<newSafety)
      {
        newSafety = distance;
      }
    }
  }
  return newSafety;
} 
