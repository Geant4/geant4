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
//
//
//---------------------------------------------------------------
//
// G4VAuxiliaryTrackInformation
//
// Class Description:
//
//  Abstract class which the user can derive his/her own concrete
// class for toring user's information associating with a G4Track
// class object.
//
//  It is user's responsibility 
//   1) Construct a concrete class object and set the pointer to
//     proper G4Track object
//   2) Concrete class derived from this class is expected to use G4Allocator
//     for memory management or something equivarent for performance reason
//
//  To set a pointer of a concrete class object to G4Track in
// G4AuxiliaryTrackingInformation concrete implementation, given the G4Track
// object is available only by "pointer to const", SetAuxiliaryTrackInformation()
// method of G4TrackingManager is available.
//
//  The concrete class object is deleted by Geant4 kernel when
// associated G4Track object is deleted.


#ifndef G4VAuxiliaryTrackInformation_H
#define G4VAuxiliaryTrackInformation_H 1

#include "globals.hh"

class G4VAuxiliaryTrackInformation
{ 
  public: // With Description
    G4VAuxiliaryTrackInformation();
    virtual ~G4VAuxiliaryTrackInformation();

    virtual void Print() const;
};

#endif

