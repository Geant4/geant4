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
// $Id: G4HadronValues.hh,v 1.9 2005/11/23 11:24:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

//
//  G4HadronValues header file
//
//
//  Kinematic and dynamic values 
//  N.  Starkov 2003.
//
//  Modifications:
//  14.11.05 Use PDG code instead of static particle pointers (N.Starkov)
//  23.11.05 cleanup (V.Ivanchenko)
//

#ifndef  G4HadronValues_h
#define  G4HadronValues_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"

#define MyPi      3.141593
#define MbToGeV2  2.568
#define GeV2ToMb  0.38939
#define MbToFm2   25.68

class G4HadronValues 
{
public:
      
  G4HadronValues();
  ~G4HadronValues(); 

  void GetHadronValues(const G4DynamicParticle * aHadron);
  
  G4double  HadrTot, HadrSlope, HadrReIm,  DDSect2, DDSect3,
            MomentumCM;

};

#endif

 
