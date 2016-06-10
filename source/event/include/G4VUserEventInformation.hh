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
// $Id: G4VUserEventInformation.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
//---------------------------------------------------------------
//
// G4VUserEventInformation
//
// Class Description:
//
//  Abstract class which the user can derive his/her own concrete
// class for toring user's information associating with a G4Event
// class object.
//
//  It is user's responsibility to construct a concrete class object 
// and set the pointer to proper G4Event object.
//
//  To set a pointer of a concrete class object to G4Event in
// G4UserEventingAction concrete implementation, given the G4Event
// object is available only by "pointer to const", SetUserEventInformation()
// method of G4EventManager is available.
//  Alternatively, the user may modify GenerateEvent() method of 
// his/her own RunManager. 
//
//  The concrete class object is deleted by Geant4 kernel when
// associated G4Event object is deleted.


#ifndef G4VUserEventInformation_H
#define G4VUserEventInformation_H 1

class G4VUserEventInformation
{
  public:
    G4VUserEventInformation() {;}
    virtual ~G4VUserEventInformation() {;}

  public:
    virtual void Print() const = 0;
};

#endif

