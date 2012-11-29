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
// $Id$
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestStreamLogger
//
// Class description:
//
// A G4GeomTestLogger that outputs to a stream.

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestStreamLogger_hh
#define G4GeomTestStreamLogger_hh

#include "G4Types.hh"
#include "G4GeomTestLogger.hh"

class G4VPhysicalVolume;

class G4GeomTestStreamLogger : public G4GeomTestLogger
{
  public:  // with description

    G4GeomTestStreamLogger( std::ostream &o,
                            G4int theMaxPointsPerError=20 );
    virtual ~G4GeomTestStreamLogger();
      // Constructors and virtual destructor

    virtual void SolidProblem( const G4VSolid *solid, 
                               const G4String &message,
                               const G4ThreeVector &point );
    virtual void OverlappingDaughters( const G4GeomTestOverlapList *list );
    virtual void OvershootingDaughter( const G4GeomTestOvershootList *list );
    virtual void NoProblem ( const G4String &message );

  protected:
  
    void PrintSegmentListHeader();
    void PrintSegmentListElement( const G4ThreeVector &s1,
                                  const G4ThreeVector &s2 );
    // Format aides

    const char *IsAre(G4int n);
  
  public:  // Solaris 7 insists on these guys being public

    class PrintPos
    {
      public:
        PrintPos(const G4ThreeVector pos, G4bool useUnit=true )
          : p(pos), unit(useUnit) {;}
    
        void Print( std::ostream & ) const;
  
      private:
        G4ThreeVector p;
        G4bool unit;
    };
  
    class VolumeNameAndCopy
    {
      public:
        VolumeNameAndCopy( const G4VPhysicalVolume *theVolume )
          : volume(theVolume) {;}
    
        void Print( std::ostream & ) const;
  
      private:
        const G4VPhysicalVolume *volume;
    };
  
    friend std::ostream &operator<<( std::ostream &,
                           const G4GeomTestStreamLogger::PrintPos & );
    friend std::ostream &operator<<( std::ostream &,
                           const G4GeomTestStreamLogger::VolumeNameAndCopy & );
  
  protected:

    std::ostream &out;
    G4int maxPointsPerError;
};

#endif
