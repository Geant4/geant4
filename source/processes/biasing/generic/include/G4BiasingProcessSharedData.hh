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
// $Id: $
//
//--------------------------------------------------------------------------
//
// G4BiasingProcessSharedData
//
// Class Description:
//        This class represents the data that the G4BiasingProcessInterface
//    objects attached to a same particle / hold by a same G4ProcessManager
//    share. It allows the active G4BiasingProcessInterface objects to share
//    information on operations that are common the these instances:
//        - the current and previous G4VBiasingOperator objects (their
//          pointers are collected once, at the beginning of the PostStepGPIL
//          by the first interface invoked)
//        - add the list of cooperating G4BiasingProcessInterface objects,
//          acting on a same particle type.
//        G4BiasingProcessInterface is friend class of this one.
//---------------------------------------------------------------------------
//   Initial version                         Sep. 2014 M. Verderi

#ifndef G4BiasingProcessSharedData_h
#define G4BiasingProcessSharedData_h

#include "globals.hh"
#include "G4Cache.hh"
#include <vector>

class G4VBiasingOperator;
class G4BiasingProcessInterface;
class G4ProcessManager;
class G4ParallelGeometriesLimiterProcess;

class G4BiasingProcessSharedData {

friend class G4BiasingProcessInterface;
friend class G4ParallelGeometriesLimiterProcess;

public:
  // -------------------------
  // -- Public access methods:
  // -------------------------
  // -- The biasing process interface objects sharing this shared data class:
  const std::vector< const G4BiasingProcessInterface* >&           GetBiasingProcessInterfaces() const
  { return           fPublicBiasingProcessInterfaces; }
  const std::vector< const G4BiasingProcessInterface* >&    GetPhysicsBiasingProcessInterfaces() const
  { return    fPublicPhysicsBiasingProcessInterfaces; }
  const std::vector< const G4BiasingProcessInterface* >& GetNonPhysicsBiasingProcessInterfaces() const
  { return fPublicNonPhysicsBiasingProcessInterfaces; }
  
  // -- The possible geometry limiter process:
  const G4ParallelGeometriesLimiterProcess*                GetParallelGeometriesLimiterProcess() const
  { return fParallelGeometriesLimiterProcess; }
  

private:
  // -- Methods used by the G4BiasingProcessInterface objects, thanks to class friendness.
  // -- Object is created by G4BiasingProcessInterface object:
  G4BiasingProcessSharedData( const G4ProcessManager* mgr)
    : fProcessManager                  (mgr),
      fCurrentBiasingOperator          ( nullptr ),
      fPreviousBiasingOperator         ( nullptr ),
      fParallelGeometryOperator        ( nullptr ),
      fMassGeometryOperator            ( nullptr ),
      fIsNewOperator                   (true),
      fLeavingPreviousOperator         (false),
      fParallelGeometriesLimiterProcess( nullptr )
  {}
  ~G4BiasingProcessSharedData() {}
  // -- biasing operators:
  void                 CurrentBiasingOperator( G4VBiasingOperator* );
  G4VBiasingOperator*  CurrentBiasingOperator() const;
  void                PreviousBiasingOperator( G4VBiasingOperator* );
  G4VBiasingOperator* PreviousBiasingOperator() const;

private:
  // -- 
  const G4ProcessManager* fProcessManager;
  // -- biasing operators:
  G4VBiasingOperator*   fCurrentBiasingOperator;
  G4VBiasingOperator*  fPreviousBiasingOperator;
  G4VBiasingOperator* fParallelGeometryOperator;
  G4VBiasingOperator*     fMassGeometryOperator;
  // -- 
  G4bool                         fIsNewOperator;
  G4bool               fLeavingPreviousOperator;

  // -- biasing process interfaces sharing this object:
  std::vector < G4BiasingProcessInterface* >                       fBiasingProcessInterfaces;
  std::vector < G4BiasingProcessInterface* >                fPhysicsBiasingProcessInterfaces;
  std::vector < G4BiasingProcessInterface* >             fNonPhysicsBiasingProcessInterfaces;
  // -- the same ones, for public use:
  std::vector < const G4BiasingProcessInterface* >           fPublicBiasingProcessInterfaces;
  std::vector < const G4BiasingProcessInterface* >    fPublicPhysicsBiasingProcessInterfaces;
  std::vector < const G4BiasingProcessInterface* > fPublicNonPhysicsBiasingProcessInterfaces;

  // -- possible process limiting step on parallel geometries:
  G4ParallelGeometriesLimiterProcess*                      fParallelGeometriesLimiterProcess;

  
  // -- thread local:
  // -- Map between process managers and shared data. This map is made of
  // -- pointers of G4BiasingSharedData instead of objects themselves :
  // -- each process needs to keep a valid pointer of a shared data object
  // -- but a map of object will make pointers invalid when map is increased.
  static G4MapCache< const G4ProcessManager*, 
		     G4BiasingProcessSharedData* >         fSharedDataMap;

};

#endif
