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
// Description:
//      Manage the per-thread state of solids - those which 
//        have a per-thread state and dependent classes (if any)
//   In particular it 
//       - owns the arrays that implement 'split' classes 
//       - classes/objects which are owned by the split classes.
//   Background: the classes/objects affected are  
//       - 'split' classes part of its state is per-thread,
//       - per-thread objects, in particular those which are owned 
//         by the split classes.
// Goal: Take ownership and control of per-thread state of 
//        classes to work with multi-threading. 
//       Offshoot of G4GeometryWorkspace, to deal with Particles.
// 
// Designed / created by John Apostolakis
// Interface design - review with Andrea Dotti.
// 
// First version: 4th Oct 2013
//    Created due to dependency issue with G4GeometryWorkspace 
// Working version:  

#ifndef G4PARTICLESWORKSPACE_HH
#define G4PARTICLESWORKSPACE_HH

#include "G4ParticleDefinition.hh"

#include "G4TWorkspacePool.hh"

class G4ParticlesWorkspace
{
  public:

    typedef G4TWorkspacePool<G4ParticlesWorkspace> pool_type;
      G4ParticlesWorkspace(G4bool verbose=false);
     ~G4ParticlesWorkspace();

     void UseWorkspace();     //Take ownership
     void ReleaseWorkspace(); //Release ownership
     void DestroyWorkspace(); //Release ownership and destroy

     void InitialiseWorkspace();
      // To be called at start of each run (especially 2nd and further runs)

     void   SetVerbose(G4bool v) { fVerbose=v; } 
     G4bool GetVerbose()  { return fVerbose;   } 
  
    static pool_type* GetPool();

 protected:  // Implementation methods

      void   InitialiseParticles();

 private:

     // Helper pointer - can be per instance or shared
     G4PDefManager *fpPDefSIM;
  
  // Per Instance variables
  //   NOTE: the ownership of the Data Arrays is IN this object

     // Store SubInstanceManager object pointers (SIM pointers)
     G4PDefData  *fpPDefOffset;
     G4bool       fVerbose;
};

#endif //G4PARTICLESWORKSPACE_HH
