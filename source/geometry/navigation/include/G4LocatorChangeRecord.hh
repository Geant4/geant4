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
// G4LocatorChangeRecord 
//
// Class description: 
//
// Record the changes in an endpoint of a locator.
// Its key use is in playing these back in case of a problem.

// Author: John Apostolakis, 27.08.19 - First version
// --------------------------------------------------------------------
#ifndef G4LOCATOR_CHANGE_RECORD_HH
#define G4LOCATOR_CHANGE_RECORD_HH

#include <vector>
#include "G4FieldTrack.hh"

class G4LocatorChangeRecord
{
  public:
 
    enum EChangeLocation { kInvalidCL = 0, kUnknownCL = 1,                // 2
                           kInitialisingCL, kIntersectsAF, kIntersectsFB, // 3
                           kNoIntersectAForFB,  kRecalculatedB,           // 2
                           kInsertingMidPoint,  kRecalculatedBagn,        // 2
                           kLevelPop };  
 
    static const char* fNameChangeLocation[];
    static const char* GetNameChangeLocation( EChangeLocation );
   
    G4LocatorChangeRecord( EChangeLocation codeLocation,
                           G4int iter,
                           unsigned int count,
                           const G4FieldTrack& fieldTrack )
        : fCodeLocation( codeLocation), fIteration(iter), fEventCount(count),
          fFieldTrack( fieldTrack ) {}
    G4LocatorChangeRecord( const G4LocatorChangeRecord &  ) = default;
    G4LocatorChangeRecord(       G4LocatorChangeRecord && ) = default;

    // No set methods -> create a new record for each entry (more reliable)
    // void SetLocation( EChangeLocation loc ) { fCodeLocation= loc; }
    // void SetLength( double len ) { fLength= len; }
    // void SetCount( int cnt ) {  fEventCount= cnt; }
    // void SetIteration( int iter ) {  fIteration= iter; }

    inline EChangeLocation GetLocation() const { return fCodeLocation; }
    inline unsigned int   GetCount()     const { return fEventCount; }
    inline G4int          GetIteration() const { return fIteration; }
    inline G4double GetLength() const { return fFieldTrack.GetCurveLength(); }
   
    friend std::ostream& operator<< ( std::ostream& os,
                         const G4LocatorChangeRecord& r );
      // Streaming operator, using StreamInfo().

    friend std::ostream& operator<< ( std::ostream& os,
                         const std::vector<G4LocatorChangeRecord> & vecR );
   
    std::ostream& StreamInfo(std::ostream& os) const;

    static std::ostream& ReportVector ( std::ostream& os,
                         const std::string & nameOfRecord, 
                         const std::vector<G4LocatorChangeRecord> & lcr );
   
    static std::ostream& ReportEndChanges ( std::ostream& os,
                         const std::vector<G4LocatorChangeRecord> & startA,
                         const std::vector<G4LocatorChangeRecord> & endB );
   
  private:

    EChangeLocation fCodeLocation = kInvalidCL;
    G4int        fIteration  = -1;
    unsigned int fEventCount = 0;
    G4FieldTrack fFieldTrack;
};

#endif
