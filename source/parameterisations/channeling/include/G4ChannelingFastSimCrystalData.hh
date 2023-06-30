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

#ifndef G4ChannelingFastSimCrystalData_h
#define G4ChannelingFastSimCrystalData_h 1

#include "G4ios.hh"
#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include <unordered_map>

#include "G4VChannelingFastSimCrystalData.hh"
#include "G4ChannelingFastSimInterpolation.hh"

/** \file G4ChannelingFastSimCrystalData.hh
* \brief Definition of the G4ChannelingFastSimCrystalData class
* The class inherits G4VChannelingFastSimCrystalData containing the data and properties
* related to the crystal lattice.
* The functions related to the crystal geometry (transformation of coordinates and angles
* from the reference system of the bounding box of the local volume to
* the crystal lattice co-rotating reference system and vice versa) and
* initialization function SetMaterialProperties are compiled in this class.
*/

class G4ChannelingFastSimCrystalData  : public G4VChannelingFastSimCrystalData
{
public:

    G4ChannelingFastSimCrystalData();
    virtual ~G4ChannelingFastSimCrystalData() = default;

public:

    ///find and upload crystal lattice input files, calculate all the basic values
    ///(to do only once)
    void SetMaterialProperties(const G4Material* crystal, const G4String &lattice);

    ///calculate the coordinates in the co-rotating reference system
    ///within a channel (periodic cell)
    ///(connected with crystal planes/axes either bent or straight)
    G4ThreeVector CoordinatesFromBoxToLattice(const G4ThreeVector &pos0);

    ///calculate the coordinates in the Box reference system
    ///(connected with the bounding box of the volume)
    G4ThreeVector CoordinatesFromLatticeToBox(const G4ThreeVector &pos);

    ///change the channel if necessary, recalculate x o y
    G4ThreeVector ChannelChange(G4double &x, G4double &y, G4double &z);

    ///calculate the horizontal angle in the co-rotating reference system
    ///within a channel (periodic cell)
    ///(connected with crystal planes/axes either bent or straight)
    G4double AngleXFromBoxToLattice(G4double tx, G4double z){return tx-AngleXShift(z);}

    ///calculate the horizontal angle in the Box reference system
    ///(connected with the bounding box of the volume)
    G4double AngleXFromLatticeToBox(G4double tx, G4double z){return tx+AngleXShift(z);}

    ///auxialiary function to transform the horizontal angle
    G4double AngleXShift(G4double z){return fMiscutAngle + z*fCurv;}

private:

    ///variables
    G4int fNsteps=353;//number of steps per channeling oscillation
    G4double fR0=1.1*CLHEP::fermi;//*A^(1/3) - radius of nucleus

    ///Values related to coordinate transformation
    long long int fNChannelx=0;//horizontal number of channel
                             //(either straight of bent) inside the box
    long long int fNChannely=0;//vertical number of channel (either straight of bent)
                             //inside the box; =0 in the case of planes

    ///values related to the crystal lattice
    G4int fNpointsx=0,fNpointsy=0;// number of horizontal and vertical nodes of interpolation
    G4double fDx=0, fDy=0;// channel (periodic cell)
                      //horizontal and vertical dimensions

    ///fundamental constants of material
    std::vector <G4double> fN0; // nuclear concentration
    std::vector <G4double> fU1; // amplitude of thermal oscillations
    std::vector <G4double> fZ1;//atomic number of each element
    std::vector <G4double> fAN; //atomic mass of each element
};

#endif
