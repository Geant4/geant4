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
// G4MultiLevelLocator 
//
// Class description:
// 
// Implementing the calculation of the intersection point with a boundary when
// PropagationInField is used. Derived from method LocateIntersectionPoint()
// from G4PropagatorInField, it is based on a linear method for finding the
// intersection point by means of a 'depth' algorithm in case of slow progress
// (intersection is not found after 100 trials).

// Author: Tatiana Nikitina (CERN), 27 October 2008
// ---------------------------------------------------------------------------
#ifndef G4MULTILEVELLOCATOR_HH
#define G4MULTILEVELLOCATOR_HH 1

#include "G4VIntersectionLocator.hh"

/**
 * @brief G4MultiLevelLocator implements the calculation of the intersection
 * point with a boundary when G4PropagationInField is used. Derived from method
 * LocateIntersectionPoint() from G4PropagatorInField, it is based on a linear
 * method for finding the intersection point by means of a 'depth' algorithm
 * in case of slow progress (intersection is not found after 100 trials).
 */

class G4MultiLevelLocator : public G4VIntersectionLocator
{
  public:
 
    /**
     * Constructor and Destructor.
     */
    G4MultiLevelLocator(G4Navigator *theNavigator);
    ~G4MultiLevelLocator() override;
     
    /**
     * If such an intersection exists, this method calculates the intersection
     * point of the true path of the particle with the surface of the current
     * volume (or of one of its daughters). 
     * Should use lateral displacement as measure of convergence.
     *  @param[in] curveStartPointTangent Start point tangent track.
     *  @param[in] curveEndPointTangent End point tangent track.
     *  @param[in] trialPoint Trial point.
     *  @param[out] intersectPointTangent Intersection point tangent track.
     *  @param[out] recalculatedEndPoint Flagging if end point was recomputed.
     *  @param[in,out] fPreviousSafety Previous safety distance.
     *  @param[in,out] fPreviousSftOrigin Previous safety point origin.
     *  @returns Whether intersection exists or not. 
     */
    G4bool EstimateIntersectionPoint( 
        const  G4FieldTrack&       curveStartPointTangent,           // A
        const  G4FieldTrack&       curveEndPointTangent,             // B
        const  G4ThreeVector&      trialPoint,                       // E
               G4FieldTrack&       intersectPointTangent,            // Output
               G4bool&             recalculatedEndPoint,             // Out
               G4double&           fPreviousSafety,                  // In/Out
               G4ThreeVector&      fPreviousSftOrigin) override;     // In/Out

    /**
     * Dumps statistics.
     */
    void ReportStatistics(); 

    /**
     * Setters.
     */
    inline void SetMaxSteps(unsigned int valMax) { fMaxSteps = valMax; }
    inline void SetWarnSteps(unsigned int valWarn) { fWarnSteps = valWarn; }

  private:

    void ReportFieldValue( const G4FieldTrack& locationPV,
                           const char* nameLoc,
                           const G4EquationOfMotion* equation );

    // Invariants -- parameters
    // ====================================   
    static const G4int max_depth = 10;
    unsigned int fMaxSteps = 10000; // Effort abandoned; signal is looping 
    unsigned int fWarnSteps = 1000; // Warn about many steps (but succeeded)

    // State - varies during simulation
    // ====================================
    G4FieldTrack* ptrInterMedFT[max_depth+1];  // Used to store intermediate
                                               // tracks values in case of too
                                               // slow progress
    // Counters for statistics & debugging
    // ====================================
    unsigned long int fNumCalls = 0; 
    unsigned long int fNumAdvanceFull = 0,
                      fNumAdvanceGood = 0,
                      fNumAdvanceTrials = 0;
};

#endif
