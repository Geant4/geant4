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
//    **************************************
//    *                                    *
//    *     CellSurvivalV79LQCeraD1.cc     *
//    *                                    *
//    **************************************
//
// Author: Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 12 October 2006 B. Mascialino      first implementation
// -------------------------------------------------------------------


#include "globals.hh"
#include "CellSurvivalV79LQCeraD1.hh"
#include "CellPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"

 
CellSurvivalV79LQCeraD1::CellSurvivalV79LQCeraD1()
{

}

CellSurvivalV79LQCeraD1::~CellSurvivalV79LQCeraD1()
{

}

void CellSurvivalV79LQCeraD1::SurvivalFormula(G4double dose)
{
  G4double alpha = 4.00 * ( gray );
  G4double alpha_gray = alpha * gray; 
  G4double beta = 25.00 * gray;

  // DEBUGGING
  // G4cout << "I WILL SURVIVE !!!!!" << G4endl;

  //
  // ******************************************************************************
  // BIBLIOGRAPHY: Cera et al, Int J Rad Biol, vol. 69, no. 6, pp. 729-738, 1996.
  // Cell line: V79-753B
  // Particle: deuteron
  // Energy: 3.20 MeV/n
  // Survival model: Linear-quadratic
  // Parameters:
  //            alpha = (0.25 +/- 0.04) Gy^-1
  //            beta  = (0.04 +/- 0.01) Gy^-2
  //********************************************************************************
  //

  probability = exp ( - ( ( dose/gray ) / ( alpha_gray / ( gray * gray ) ) + ( dose/gray ) / ( beta / gray ) ) );


  // DEBUGGING
  // G4cout << "alpha_gray= " <<  1 / ( alpha_gray / ( gray * gray ) )  << G4endl;
  // G4cout << "beta = " <<   1 / ( beta / gray )  << G4endl;
  // G4cout << "alpha_gray * dose/gray * dose/gray= " <<  ( ( dose * dose ) / ( gray * gray ) ) / ( alpha_gray / ( gray * gray ) ) << G4endl;
  // G4cout << "beta * dose/gray= " <<  ( dose/gray ) / ( beta / gray ) << G4endl;
  // G4cout << "alpha_gray * dose/gray * dose/gray + beta * dose/gray= " << ( dose/gray ) / ( alpha_gray / ( gray * gray ) ) + ( dose/gray ) / ( beta / gray )   << G4endl;
  // G4cout << "probability ( debugging ) = " <<  exp ( - ( ( dose/gray ) / ( alpha_gray / ( gray * gray ) ) + ( dose/gray ) / ( beta / gray ) ) )   << G4endl;



 CellPrimaryGeneratorAction* primary = 
   (CellPrimaryGeneratorAction*) G4RunManager::GetRunManager()-> GetUserPrimaryGeneratorAction();

  // Type of the primary particle
 G4String primaryParticleName = primary -> GetParticle();
   
  // Energy of the primary particle
 G4double primaryParticleEnergy = primary -> GetInitialEnergy();

  // TESTING
  // Expected survival probability
  // G4cout << "Expected survival probability= " << 0.99999 << G4endl;
 

 G4cout << "Primary particle: " << primaryParticleName << " with energy " <<
   primaryParticleEnergy/MeV << " MeV" << G4endl;
 

}
G4double CellSurvivalV79LQCeraD1::GetSurvival()
{
  return probability; 
}
