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
// $Id: G4VUserRegionInformation.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//
// G4VUserRegionInformation
//
// Class Description:
//
// Abstract class, the user can subclass for storing information
// associated with a G4Region class object.
//
// It is user's responsibility to construct a concrete class object 
// and set the pointer to proper G4Region object.
// The concrete class object is deleted by Geant4 kernel when
// associated G4Region object is deleted.

// History:
// 21/10/2003 - M.Asai, Created.
// --------------------------------------------------------------------
#ifndef G4VUserRegionInformation_H
#define G4VUserRegionInformation_H 1

class G4VUserRegionInformation
{
  public:  // with description
  
    G4VUserRegionInformation() {}
    virtual ~G4VUserRegionInformation() {}

    virtual void Print() const = 0;
};

#endif

