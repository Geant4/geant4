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
// G4PhysicsVectorType
//
// Description:
//
// Enumerator to define the physics vector type:
//   G4PhysicsFreeVector
//   G4PhysicsLinearVector
//   G4PhysicsLogVector
//
// Spline type:
//   G4SplineType::Simple     - 2d derivative continues
//   G4SplineType::Base       - 3d derivative continues
//   G4SplineType::FixedEdges - 3d derivatives continues, 1st and last 
//                              derivatives are fixed 
//
// Author: H.Kurashige, 9 March 2001
//
// --------------------------------------------------------------------
#ifndef G4PhysicsVectorType_hh
#define G4PhysicsVectorType_hh 1

enum G4PhysicsVectorType
{
  T_G4PhysicsFreeVector = 0,
  T_G4PhysicsLinearVector,
  T_G4PhysicsLogVector
};

enum class G4SplineType
{
  Simple = 0,
  Base,
  FixedEdges
};

#endif
