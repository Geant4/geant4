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
// $Id: G4HarmonicPolMagField.hh 68055 2013-03-13 14:43:28Z gcosmo $
//
// class G4HarmonicPolMagField
//
// Class description:
//
// Class describing magnetic field parametrised by harmonic polynom up to
// 3rd order. The function MagneticField(yTrack,B) calculates the magnetic
// field induction vector B for the trajectory point yTrack according to
// formula given in:
//   M.Metcalf, Analysis of the SFM Field,
//              OM Development Note AP-10 (revised), 1974

// History:
// 3.2.97 - V.Grichine, created.
// --------------------------------------------------------------------

#ifndef G4HARMONICPOLMAGFIELD_HH
#define G4HARMONICPOLMAGFIELD_HH

#include "G4MagneticField.hh"

class G4HarmonicPolMagField : public G4MagneticField
{
  public:  // with description
                       
    G4HarmonicPolMagField();
   ~G4HarmonicPolMagField();
     
    void GetFieldValue(const G4double yTrack[] ,
                             G4double B[]      ) const  ;
    G4HarmonicPolMagField* Clone() const;
};

#endif
