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
// $Id: G4VUserPrimaryParticleInformation.hh 67971 2013-03-13 10:13:24Z gcosmo $
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
    G4VUserPrimaryParticleInformation();
    virtual ~G4VUserPrimaryParticleInformation();

  public:
    virtual void Print() const = 0;
};

#endif

