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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//    Extends G4Track properties with information needed for handling
//    in biasing.
//    These information regard at this point under what conditions
//    the track was borned.
//
//      ----------------G4BiasingTrackData ----------------
//
// Author: M.Verderi (LLR), November 2013
//
// --------------------------------------------------------------------

#ifndef G4BiasingTrackData_hh
#define G4BiasingTrackData_hh

class G4Track;
class G4VBiasingOperation;
class G4VBiasingOperator;
class G4BiasingProcessInterface;

class G4BiasingTrackData {
public:
  G4BiasingTrackData(const G4Track* track);
  G4BiasingTrackData(const G4Track*                            track,
		     const G4VBiasingOperation*       birthOperation,
		     const G4VBiasingOperator*         birthOperator,
		     const G4BiasingProcessInterface*   birthProcess);
  ~G4BiasingTrackData();
  
  void SetBirthOperation( const G4VBiasingOperation*       birthOperation ) { fBirthOperation  = birthOperation; }
  void SetBirthOperator ( const G4VBiasingOperator*         birthOperator ) { fBirthOperator   = birthOperator;  }
  void SetBirthProcess  ( const G4BiasingProcessInterface*   birthProcess ) { fBirthProcess    = birthProcess;   }
  
  const G4Track*                              GetTrack() const { return fTrack; }
  const G4VBiasingOperation*         GetBirthOperation() const { return fBirthOperation; }
  const G4VBiasingOperator*           GetBirthOperator() const { return fBirthOperator; }
  const G4BiasingProcessInterface*     GetBirthProcess() const { return fBirthProcess; }
  
private:
  const G4Track*                           fTrack;
  const G4VBiasingOperation*      fBirthOperation;
  const G4VBiasingOperator*        fBirthOperator;
  const G4BiasingProcessInterface*  fBirthProcess;
		     
};

#endif
