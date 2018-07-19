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
// $Id: G4ElasticHNScattering.hh 100828 2016-11-02 15:25:59Z gcosmo $

#ifndef G4ElasticHNScattering_h
#define G4ElasticHNScattering_h 1

// ------------------------------------------------------------
//                    GEANT 4 class header file
//
//      ---------------- G4ElasticHNScattering --------------
//                   by V. Uzhinsky, March 2008.
//             elastic scattering used by Fritiof model
//                 Take a projectile and a target
//                 scatter the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
#include "G4FTFParameters.hh"
#include "G4ThreeVector.hh"
#include "G4SampleResonance.hh"   // Uzhi 2014

class G4VSplitableHadron;
class G4ExcitedString;


class G4ElasticHNScattering {
  public:
    G4ElasticHNScattering();
    virtual ~G4ElasticHNScattering();
    virtual G4bool ElasticScattering( G4VSplitableHadron* aPartner, G4VSplitableHadron* bPartner,
                                      G4FTFParameters* theParameters ) const;

  private:
    G4ElasticHNScattering( const G4ElasticHNScattering& right );
    G4ThreeVector GaussianPt( G4double AveragePt2, G4double maxPtSquare ) const;
    const G4ElasticHNScattering& operator=( const G4ElasticHNScattering& right );
    int operator==( const G4ElasticHNScattering& right ) const;
    int operator!=( const G4ElasticHNScattering& right ) const;
};

#endif

