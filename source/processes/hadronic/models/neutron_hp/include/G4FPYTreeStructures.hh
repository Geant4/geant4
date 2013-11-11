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

