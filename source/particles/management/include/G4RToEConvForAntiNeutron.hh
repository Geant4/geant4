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
// $Id: G4RToEConvForAntiNeutron.hh,v 1.2 2002-12-16 11:15:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
// Class Description
//  This class is a Range to Energy Converter for anti Neutron.
//
// ------------------------------------------------------------
//   First Implementation          5 Oct. 2002  H.Kurahige
// ------------------------------------------------------------

#ifndef G4RToEConvForAntiNeutron_h
#define G4RToEConvForAntiNeutron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

#include "G4VRangeToEnergyConverter.hh"


class G4RToEConvForAntiNeutron : public G4VRangeToEnergyConverter
{
  public: 
  //  constructor
  G4RToEConvForAntiNeutron();

  public:
  //  destructor
  virtual ~G4RToEConvForAntiNeutron();


 // calculate energy cut from given range cut for the material
  virtual G4double Convert(G4double rangeCut, const G4Material* material);

};

inline
  G4double G4RToEConvForAntiNeutron::Convert(G4double , const G4Material* )
{
  // reutrn lowest energy
  return LowestEnergy;
}
#endif









