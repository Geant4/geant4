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
// $Id: G4VMaterialMap.hh,v 1.2 2001/07/11 10:02:25 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef G4VMaterialMap_hh
#define G4VMaterialMap_hh 1

// class description:
//
//  This is an abstract base class for bookkeeping the material name
// for retrieving persistent geometry.
// User must provide the appropriate class for constructing material
// objects from material name when retrieving a stored persistent
// geometry.
//  The virtual method LookUp() will be invoked from G4PersistentGeomMan.
//

#include "globals.hh"
#include "G4Material.hh"

class G4VMaterialMap 
{
  public:
      static G4VMaterialMap* GetMaterialMap();
      //  Static method to return the pointer to the singleton object.
      // Note that this method does NOT create the singleton object.

  protected:
      G4VMaterialMap();

  public:
      virtual ~G4VMaterialMap();

  private: 
      static G4VMaterialMap * fMaterialMap;

  public: // with description
      virtual G4Material* LookUp(const G4String aName)=0;
      //  Returns a pointer to the material

};

#endif
