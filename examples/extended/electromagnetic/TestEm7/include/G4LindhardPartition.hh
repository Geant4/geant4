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
/*
 *  \file electromagnetic/TestEm7/include/G4LindhardPartition.hh
 *  \brief Definition of the G4LindhardPartition class
 *
 *  Created by Marcus Mendenhall on 1/14/08.
 *  2008 Vanderbilt University, Nashville, TN, USA.
 *
 */

//

#include "globals.hh"

class G4Material;

class G4VNIELPartition 
{
public:
        G4VNIELPartition() { }
        virtual ~G4VNIELPartition() { }
        
        // return the fraction of the specified energy which will be deposited as NIEL
        // if an incoming particle with z1, a1 is stopped in the specified material
        // a1 is in atomic mass units, energy in native G4 energy units.
        virtual G4double PartitionNIEL(
                G4int z1, G4double a1, const G4Material *material, G4double energy
        ) const =0;
};

class G4LindhardRobinsonPartition : public G4VNIELPartition
{
public:
        G4LindhardRobinsonPartition();
        virtual ~G4LindhardRobinsonPartition() { }
        
        virtual G4double PartitionNIEL(
                G4int z1, G4double a1, const G4Material *material, G4double energy
        ) const ;
        
        G4double z23[120];
        size_t   max_z;
};

