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
// $Id: G4ProcessAttribute.hh,v 1.4 2001-07-11 10:08:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4VProcess.hh"
#include "G4ProcessManager.hh"

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

    G4ProcessAttribute & operator=(G4ProcessAttribute &right);
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
	 isActive(true)
{
  idxProcessList = -1;
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





