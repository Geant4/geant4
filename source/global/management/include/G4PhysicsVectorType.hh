// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsVectorType.hh,v 1.1 2001-03-09 03:39:27 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------
//
// G4PhysicsVectorType.hh
//
// Class Description:
//   This is an enumerator to define the physics vector type
//
#ifndef G4PhysicsVectorType_h
#define G4PhysicsVectorType_h

enum G4PhysicsVectorType
{
  T_G4PhysicsVector =0,
  T_G4PhysicsLinearVector,
  T_G4PhysicsLogVector,
  T_G4PhysicsLnVector,
  T_G4PhysicsFreeVector,
  T_G4PhysicsOrderedFreeVector,
  T_G4LPhysicsFreeVector
};

#endif
