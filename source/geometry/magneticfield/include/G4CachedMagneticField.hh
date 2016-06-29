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
// $ Id: G4CachedMagneticField.hh,v 1.0 2009/07/20 18:53:00 japost Exp $
//
//
// class G4CachedMagneticField
//
// Class description:
//
// Caches Magnetic Field value, for field whose evaluation is expensive.

// History:
// - Created. JA, July 20th, 2009.
// --------------------------------------------------------------------

#ifndef G4CACHED_MAGNETIC_FIELD_DEF
#define G4CACHED_MAGNETIC_FIELD_DEF

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

class G4CachedMagneticField : public G4MagneticField
{
  public:  // with description

     G4CachedMagneticField(G4MagneticField *, G4double distanceConst);
     virtual ~G4CachedMagneticField();
       // Constructor and destructor. No actions.

     G4CachedMagneticField(const G4CachedMagneticField &r);
     G4CachedMagneticField& operator = (const G4CachedMagneticField &p);
       // Copy constructor & assignment operator.

     virtual void  GetFieldValue( const G4double Point[4],
                                        G4double *Bfield ) const;
     
     G4double GetConstDistance() const         { return fDistanceConst; } 
     void     SetConstDistance( G4double dist ){ fDistanceConst= dist;}

     G4int GetCountCalls() const { return fCountCalls; }
     G4int GetCountEvaluations() const { return fCountEvaluations; } 
     void  ClearCounts() { fCountCalls = 0; fCountEvaluations=0; }
     void  ReportStatistics();
    
    virtual G4Field* Clone() const;

  private: 
     G4MagneticField *fpMagneticField;
     // When the field is evaluated within this distance it will not change
     G4double      fDistanceConst;
     // Caching state
     mutable G4ThreeVector fLastLocation;
     mutable G4ThreeVector fLastValue;

  protected:
     mutable G4int fCountCalls, fCountEvaluations;  
};

#endif /* G4CACHED_MAGNETIC_FIELD_DEF */
