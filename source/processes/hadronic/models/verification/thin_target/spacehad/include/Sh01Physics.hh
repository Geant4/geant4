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
// 
//      GEANT 4 class 
//
//    Sh01Physics
//              
// 
//    Modified:
//
//    06.03.03 V. Grichine (based on test30)of V. Ivanchenko)
//

#ifndef Sh01Physics_h
#define Sh01Physics_h 1

#include "globals.hh"
#include "Sh01HadronProduction.hh"

class G4VProcess;
class G4Material;

////////////////////////////////////////////////////////////////////////

class Sh01Physics 
{
  public:
  
    Sh01Physics();
   ~Sh01Physics();
     
    G4VProcess* GetProcess(const G4String&, const G4String&, G4Material*);     
    G4double GetNucleusMass() {return theProcess->GetMass();};
		                      
  private:

    void Initialise();
	     
    Sh01HadronProduction* theProcess; 
	
};

#endif

 


