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
// $Id: G4VUserRegionInformation.hh,v 1.1 2003-10-21 20:22:06 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4VUserRegionInformation
//
// Class Description:
//
//  Abstract class which the user can derive his/her own concrete
// class for toring user's information associating with a G4Region
// class object.
//
//  It is user's responsibility to construct a concrete class object 
// and set the pointer to proper G4Region object.
//
//  The concrete class object is deleted by Geant4 kernel when
// associated G4Region object is deleted.


#ifndef G4VUserRegionInformation_H
#define G4VUserRegionInformation_H 1

class G4VUserRegionInformation
{
  public:
    G4VUserRegionInformation() {;}
    virtual ~G4VUserRegionInformation() {;}

  public:
    virtual void Print() const = 0;
};

#endif

