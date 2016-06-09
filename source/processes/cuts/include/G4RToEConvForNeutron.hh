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
// $Id: G4RToEConvForNeutron.hh,v 1.2 2006/06/29 19:30:00 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
// Class Description
//  This class is a Range to Energy Converter for Neutron.
//
// ------------------------------------------------------------
//   First Implementation          5 Oct. 2002  H.Kurahige
// ------------------------------------------------------------

#ifndef G4RToEConvForNeutron_h
#define G4RToEConvForNeutron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VRangeToEnergyConverter.hh"


class G4RToEConvForNeutron : public G4VRangeToEnergyConverter
{
  public: 
  //  constructor
  G4RToEConvForNeutron();

  public:
  //  destructor
  virtual ~G4RToEConvForNeutron();

 // calculate energy cut from given range cut for the material
  virtual G4double Convert(G4double rangeCut, const G4Material* material);

};

inline
  G4double G4RToEConvForNeutron::Convert(G4double , const G4Material* )
{
  // reutrn lowest energy
  return LowestEnergy;
}
#endif









