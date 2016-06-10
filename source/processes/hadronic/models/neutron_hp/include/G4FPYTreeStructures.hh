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
 * File:   G4FPYTreeStructures.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on September 8, 2011, 10:43 AM
 */

#ifndef G4FPYTREESTRUCTURES_HH
#define	G4FPYTREESTRUCTURES_HH

#include "G4Ions.hh"
#include "globals.hh"

/** ProbabilityBranch is a tree hierarchy storage method. Each 'branch'
 *  stores information for the isotope/isomer that is stores, as well as the
 *  probability segment in the (0,1] range references it for use in random
 *  fission product sampling.
 *  \n ProbabilityBranch branches 'up' and 'down' to other branches that
 *  have a probability segment either higher or lower in the (0,1] range.
 */
struct ProbabilityBranch
{
    /** Pointer to the \p G4Ions definition that contains the data for
     *  the isotope/isomer for this ProbabilityBranch
     */
    G4Ions* Particle;

    /** Number of different incident energies are available */
    G4int IncidentEnergiesCount;
    /** The different energies available */
    G4double* IncidentEnergies;
    /** The upper limit of the probability segment, indexed by incident energy*/
    G4double* ProbabilityRangeBottom;
    /** The lower limit of the probability segment, indexed by incident energy */
    G4double* ProbabilityRangeTop;

    /** Pointer to the branch that has a higher probability */
    ProbabilityBranch* Right;
    /** Pointer to the branch that has a lower probability */
    ProbabilityBranch* Left;
};

/** ProbabilityTree is the 'root' for each tree hierarchy, and contains
 *  pointer to the first ProbabilityBranch, or 'trunk', or each tree.
 *  \n ProbabilityTree also stores the value of the highest probability
 *  segment that is stored in the tree. This enables trees that are out of
 *  range to be easily skipped over.
 */
struct ProbabilityTree
{
    /** Number of different incident energies are available */
    static G4int IncidentEnergiesCount;
    /** The different energies available */
    static G4double* IncidentEnergies;

    /** First branch, or 'trunk', in the tree */
    ProbabilityBranch* Trunk;
    /** The value of the highest probability segment stored in this tree */
    G4double* ProbabilityRangeEnd;
    /** Counter for the number of branches that this tree has **/
    G4int BranchCount;
    /** Variable to identify if this is the last tree in or not */
    G4bool IsEnd;
};

#endif	/* G4FPYTREESTRUCTURES_HH */

