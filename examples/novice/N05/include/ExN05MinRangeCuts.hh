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
// $Id: ExN05MinRangeCuts.hh,v 1.6 2006/06/29 17:52:41 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//                  14 Aug. 1998  H.Kurashige
// ------------------------------------------------------------

#ifndef ExN05MinRangeCuts_h
#define ExN05MinRangeCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "ExN05SpecialCuts.hh"


class ExN05MinRangeCuts : public ExN05SpecialCuts
{
  public:     

     ExN05MinRangeCuts(const G4String& processName ="ExN05MinRangeCuts" );

     virtual ~ExN05MinRangeCuts();

     // PostStep GPIL
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );
            
			    
  private:
  
  // hide assignment operator as private 
      ExN05MinRangeCuts(ExN05MinRangeCuts&);
      ExN05MinRangeCuts& operator=(const ExN05MinRangeCuts& right);

};

#endif

