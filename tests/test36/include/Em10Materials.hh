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
//      ---------- Em10Materials-------
//    Originally Created in Test30 by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified for TestEm10 by V. Grichine, 30 Jan 2006 
//
//

#ifndef Em10Materials_h
#define Em10Materials_h 1

#include "globals.hh"

class G4Material;

class Em10Materials
{
  public:
  
    Em10Materials();
   ~Em10Materials();
     
	G4Material* GetMaterial(const G4String&);     
                      
  private:

	void Initialise();
	     
};

#endif

 


