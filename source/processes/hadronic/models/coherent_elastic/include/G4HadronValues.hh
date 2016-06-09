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
// $Id: G4HadronValues.hh,v 1.10 2006/06/29 20:09:05 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

 
