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
//     Sh01Material 
//             
// 
//    Modified:
//   
//    06.03.03 V. Grichine (based on test30 of V. Ivanchenko)
//

#ifndef Sh01Material_h
#define Sh01Material_h 1

#include "globals.hh"

class G4Material;

////////////////////////////////////////////////////////////////////

class Sh01Material 
{
  public:
  
    Sh01Material();
   ~Sh01Material();
     
    G4Material* GetMaterial(const G4String&);     
                      
  private:

   void Initialise();
	     
};

#endif

 


