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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test19Material -------
// 
//    Converted from Test29 to Test19 by Mikhail Kossov, 29 Jan 2005 
//
//===========================================================================

#ifndef Test19Material_h
#define Test19Material_h 1

#include "globals.hh"

class G4Material;

class Test19Material 
{
  public:
  
    Test19Material();
   ~Test19Material();
     
	G4Material* GetMaterial(const G4String&);     
                      
  private:

	void Initialise();
	     
};

#endif

 


