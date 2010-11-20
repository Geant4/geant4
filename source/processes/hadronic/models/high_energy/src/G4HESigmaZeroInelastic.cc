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
// $Id: G4HESigmaZeroInelastic.cc,v 1.12 2010-11-20 04:01:33 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G4ios.hh"

// G4 Process: Gheisha High Energy Collision model.
// This includes the high energy cascading model, the two-body-resonance model
// and the low energy two-body model. Not included are the low energy stuff
// like nuclear reactions, nuclear fission without any cascading and all
// processes for particles at rest.  
// First work done by J.L.Chuma and F.W.Jones, TRIUMF, June 96.  
// H. Fesefeldt, RWTH-Aachen, 23-October-1996
// Last modified: 29-July-1998 
 
#include "G4HESigmaZeroInelastic.hh"
#include "G4Gamma.hh"

G4HadFinalState*
G4HESigmaZeroInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                      G4Nucleus& targetNucleus)
{
  G4HEVector* pv = new G4HEVector[MAXPART];
  const G4HadProjectile* aParticle = &aTrack;
  G4HEVector incidentParticle(aParticle);
     
  G4HELambdaInelastic theLambdaInelastic;
  theLambdaInelastic.SetMaxNumberOfSecondaries(MAXPART);
  theLambdaInelastic.SetVerboseLevel(verboseLevel);
    
  G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
  G4double pgam = G4UniformRand()*incidentTotalMomentum*0.75;
  G4HEVector incidentLambda; 
  incidentLambda.SmulAndUpdate(incidentParticle, 
                      (incidentTotalMomentum - pgam)/incidentTotalMomentum);                    
  G4DynamicParticle* aLambda = new G4DynamicParticle();
  aLambda->SetDefinition(G4Lambda::Lambda());
  aLambda->SetMomentum(incidentLambda.getMomentum());
  G4HadProjectile aLambdaTrack(*aLambda);
  G4HadFinalState* result =
          theLambdaInelastic.ApplyYourself(aLambdaTrack, targetNucleus);         	
  vecLength = theLambdaInelastic.GetNumberOfSecondaries();
    
  pv[vecLength] = Gamma;
  pv[vecLength].setMomentum(incidentParticle.getMomentum());
  pv[vecLength].SmulAndUpdate( pv[vecLength],pgam/incidentTotalMomentum);
  G4DynamicParticle* aPhoton = new G4DynamicParticle();
  aPhoton->SetDefinition(G4Gamma::Gamma());
  aPhoton->SetMomentum(pv[vecLength].getMomentum());
  result->AddSecondary(aPhoton);
  delete [] pv;
  return result;
} 

