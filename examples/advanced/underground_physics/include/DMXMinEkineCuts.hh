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
// $Id: DMXMinEkineCuts.hh,v 1.2 2002-06-18 10:17:30 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//                  14 Aug. 1998  H.Kurashige
// ------------------------------------------------------------

#ifndef DMXMinEkineCuts_h
#define DMXMinEkineCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "DMXSpecialCuts.hh"


class DMXMinEkineCuts : public DMXSpecialCuts
{
  public:     

     DMXMinEkineCuts(const G4String& processName ="DMXMinEkineCuts" );

     virtual ~DMXMinEkineCuts();

     // PostStep GPIL
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );
            
			    
  private:
  
  // hide assignment operator as private 
      DMXMinEkineCuts(DMXMinEkineCuts&);
      DMXMinEkineCuts& operator=(const DMXMinEkineCuts& right);

};

#endif

