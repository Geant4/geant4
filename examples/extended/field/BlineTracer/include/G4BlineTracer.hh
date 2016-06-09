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
// $Id: G4BlineTracer.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 
// --------------------------------------------------------------------
//
// G4BlineTracer
//
// Class description:
//
// Defines a tool to trace and visualise magnetic field lines
// To use this tool in a Geant4 application the user should  
// create an instance of this class in the code as a run action.
// It will only work if a G4MagneticField field object is declared.

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlineTracer_h
#define G4BlineTracer_h 1

#include <vector>

#include "G4Types.hh"
#include "G4UserRunAction.hh"

class G4VUserPrimaryGeneratorAction;
class G4MagneticField;
class G4FieldManager;
class G4ChordFinder;

class G4BlineTracerMessenger;
class G4BlineSteppingAction;
class G4BlineEventAction;
class G4BlinePrimaryGeneratorAction;
class G4BlineEquation;

class G4BlineTracer : public G4UserRunAction 
{
   public:  // with description
  
     G4BlineTracer();
     virtual ~G4BlineTracer();
  
     virtual void BeginOfRunAction(const G4Run* aRun);
     virtual void EndOfRunAction(const G4Run* aRun);

     void ComputeBlines(G4int nlines);

     inline void SetMaxTrackingStep(G4double max_step)
       { MaxTrackingStep=max_step; }
     inline G4BlineEventAction* GetEventAction()
       { return fEventAction; }

   private:

     void ResetChordFinders();

   private:

     G4BlineTracerMessenger* fMessenger;
     G4BlineSteppingAction* fSteppingAction;
     G4BlineEventAction* fEventAction;
     G4BlinePrimaryGeneratorAction* fPrimaryGeneratorAction;
     G4double MaxTrackingStep;
     G4bool was_ResetChordFinders_already_called;
 
     G4VUserPrimaryGeneratorAction* fUserPrimaryAction;
       // User defined primary generator action

     std::vector<G4ChordFinder* > vecChordFinders;
     std::vector<G4FieldManager* > vecFieldManagers;
     std::vector<G4MagneticField* > vecMagneticFields;
     std::vector<G4BlineEquation* > vecEquationOfMotion;
       // ChordFinders, detector fields, equation of motions, and field
       // manager for the different local and global magnetic fields.
};

#endif
