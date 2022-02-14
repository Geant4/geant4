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
// G4LocatorChangeLogger
//
// Class description: 
//
// Aggregate the records of changes in an endpoint of a locator.
// Its key use is in playing these back in case of a problem.

// Author: John Apostolakis, 04.09.19 - First version
// --------------------------------------------------------------------
#ifndef G4LOCATOR_CHANGE_LOGGER_HH
#define G4LOCATOR_CHANGE_LOGGER_HH

#include <vector>
#include "G4LocatorChangeRecord.hh"
#include "G4FieldTrack.hh"

class G4LocatorChangeLogger : public std::vector<G4LocatorChangeRecord>
{
  public:

    G4LocatorChangeLogger( const std::string name ) : fName(name) {}

    void AddRecord(       G4LocatorChangeRecord && chngRecord );
    void AddRecord( const G4LocatorChangeRecord &  chngRecord );

    // Create a new record with full information
    inline
    void AddRecord( G4LocatorChangeRecord::EChangeLocation codeLocation,
                    G4int             iter,
                    unsigned int    count,
                    const G4FieldTrack &  fieldTrack );

    friend std::ostream& operator << ( std::ostream& os,
                         const G4LocatorChangeLogger& logR );
    std::ostream& StreamInfo(std::ostream& os) const;

    static std::ostream& ReportEndChanges ( std::ostream& os,
                                            const G4LocatorChangeLogger& startA,
                                            const G4LocatorChangeLogger& endB );
      // Print the changes in start, end points in columns
      // One event per row

  private:

    const std::string  fName;
};

// --------------
// Inline methods
// --------------

void G4LocatorChangeLogger::
AddRecord( G4LocatorChangeRecord::EChangeLocation codeLocation,
           G4int iter, unsigned int count,
           const G4FieldTrack &  fieldTrack )
{
  this->push_back(G4LocatorChangeRecord(codeLocation, iter, count, fieldTrack));
}

inline   
void G4LocatorChangeLogger::
AddRecord( const G4LocatorChangeRecord& chngRecord )
{
  this->push_back( chngRecord );
}

inline
void G4LocatorChangeLogger::
AddRecord( G4LocatorChangeRecord && chngRecord )
{
  this->push_back( chngRecord );
}
   
#endif
