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
// $Id: G4VUserPrimaryParticleInformation.hh,v 1.1 2003/09/12 21:51:32 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
//---------------------------------------------------------------
//
// G4VUserPrimaryParticleInformation
//
// Class Description:
//
//  Abstract class which the user can derive his/her own concrete
// class for toring user's information associating with a G4PrimaryParticle
// class object.
//
//  It is user's responsibility to construct a concrete class object 
// and set the pointer to proper G4PrimaryParticle object.
//
//  The concrete class object is deleted by Geant4 kernel when
// associated G4PrimaryParticle object is deleted.


#ifndef G4VUserPrimaryParticleInformation_H
#define G4VUserPrimaryParticleInformation_H 1

class G4VUserPrimaryParticleInformation
{
  public:
    G4VUserPrimaryParticleInformation() {;}
    virtual ~G4VUserPrimaryParticleInformation() {;}

  public:
    virtual void Print() const = 0;
};

#endif

