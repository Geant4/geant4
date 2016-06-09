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
// $Id: G4VUserEventInformation.hh,v 1.1 2003/09/09 20:09:17 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

