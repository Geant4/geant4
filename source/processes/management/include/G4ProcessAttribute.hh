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
// $Id: G4ProcessAttribute.hh 71231 2013-06-12 13:06:28Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1997, H.Kurashige
//   ----------------  G4ProcessAttribute -----------------
// Class Description
//  This class is used by G4ProcessManager ONLY for booking !!!
//
// History:
//   adds copy constructor            27 June 1998 H.Kurashige
// ------------------------------------------------------------

#ifndef G4ProcessAttribute_h
#define G4ProcessAttribute_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4ProcessManager.hh"

class G4VProcess;

class G4ProcessAttribute
{
  // this class is used by G4ProcessManager ONLY for booking !!!
  friend class G4ProcessManager;
  public:
    G4ProcessAttribute();
    G4ProcessAttribute(const G4VProcess* aProcess);
    G4ProcessAttribute(const G4ProcessAttribute &right);
    //  Constructors

    ~G4ProcessAttribute();
    //  Destructor

    G4ProcessAttribute & operator=(const G4ProcessAttribute &right);
    // Assignment operator

    G4int operator==(const G4ProcessAttribute &right) const;
    G4int operator!=(const G4ProcessAttribute &right) const;
    // equal / unequal operator

  
  protected:
    G4VProcess*           pProcess;
    // pointer to G4VProcess

    G4bool                isActive;
    // flag for activation/inactivation

    G4int                 idxProcessList;
    // index to a ProcessVector for theProcessList and 

    G4int                 idxProcVector[G4ProcessManager::SizeOfProcVectorArray];
    // index to ProcessVectors for "Doit"s and "GetPhysicalInteractionLength"s
    //   -1 : not applicable

    G4int                 ordProcVector[G4ProcessManager::SizeOfProcVectorArray];
    // ordering parameter 
};

inline 
 G4ProcessAttribute::G4ProcessAttribute(const G4VProcess* aProcess):
         pProcess((G4VProcess*)aProcess),
	 isActive(true),
	 idxProcessList(-1)
{
  for(size_t ii=0; ii<G4ProcessManager::SizeOfProcVectorArray; ii++){
    idxProcVector[ii]=-1;
    ordProcVector[ii]=0; 
  }
}

inline 
 G4int  G4ProcessAttribute::operator==(const G4ProcessAttribute &right) const
{
    return this->pProcess == right.pProcess;
}

inline 
 G4int  G4ProcessAttribute::operator!=(const G4ProcessAttribute &right) const
{
    return this->pProcess != right.pProcess;
}

#endif





