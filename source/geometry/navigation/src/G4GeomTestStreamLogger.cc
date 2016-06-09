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
// $Id: G4GeomTestStreamLogger.cc,v 1.3 2006-06-29 18:36:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestStreamLogger
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestStreamLogger.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include <iomanip>

#include "G4GeomTestOverlapList.hh"
#include "G4GeomTestOvershootList.hh"

//
// Constructor and destructor
//
G4GeomTestStreamLogger::G4GeomTestStreamLogger( std::ostream &o,
                                                G4int theMaxPointsPerError )
  : out(o), maxPointsPerError(theMaxPointsPerError)
{;}

G4GeomTestStreamLogger::~G4GeomTestStreamLogger()
{;}

//
// ::PrintPos
//
// Utility class for printing a 3 vector position
//
void G4GeomTestStreamLogger::PrintPos::Print( std::ostream &o ) const
{
  o << std::setprecision(6) << std::setw(14) << p.x()/cm;
  o << std::setprecision(6) << std::setw(14) << p.y()/cm;
  o << std::setprecision(6) << std::setw(14) << p.z()/cm;
  if (unit) o << " cm";
}

std::ostream &operator<<(      std::ostream &o,
                           const G4GeomTestStreamLogger::PrintPos &p )
{
  p.Print(o);
  return o;
}

//
// ::VolumeNameAndCopy
//
// Utility class for printing a volume's name and copy number
//
void
G4GeomTestStreamLogger::VolumeNameAndCopy::Print( std::ostream &o ) const
{
        o << volume->GetName() << "[" << volume->GetCopyNo() << "]";
}

std::ostream &operator<<( std::ostream &o,
                            const G4GeomTestStreamLogger::VolumeNameAndCopy &p )
{
  p.Print(o);
  return o;
}


//
// SolidProblem
//
void G4GeomTestStreamLogger::SolidProblem( const G4VSolid *solid, 
                                           const G4String &message,
                                           const G4ThreeVector &point )
{
  out << "GeomTest Error: SolidProblem\n"
      << "    " << message << "\n"
      << "    Solid name = " << solid->GetName() << "\n"
      << "    Local position = " << PrintPos(point) << std::endl;
}


//
// NoProblem
//
void G4GeomTestStreamLogger::NoProblem( const G4String &message )
{
  out << message << std::endl;
}


//
// OverlappingDaughters
//
void
G4GeomTestStreamLogger::OverlappingDaughters(const G4GeomTestOverlapList *list)
{
  G4int n = list->NumError();
  if (n <= 0) return;

  out << "GeomTest Error: Overlapping daughter volumes\n"
   << "    The volumes " << VolumeNameAndCopy(list->GetDaughter1())
   << " and " << VolumeNameAndCopy(list->GetDaughter2()) << ",\n"
   << "    both daughters of volume " << VolumeNameAndCopy(list->GetMother())
   << ",\n"
   << "    appear to overlap at the following " << (n>1 ? "points" : "point")
   << " in global coordinates:";

  G4int nInterval, nStop;

  if (n <= maxPointsPerError) {
    out << "\n";
    nInterval = 1;
    nStop = n;
  }
  else {
    out << " (list truncated)\n";
    nInterval = n/maxPointsPerError;
    nStop = maxPointsPerError*nInterval;
  }

  G4int i;
  G4ThreeVector s1, s2;
  
  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetGlobalPoints( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }

  out << "    Which in the mother coordinate system " << IsAre(n) << ":\n";

  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetMotherPoints( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }

  out << "    Which in the coordinate system of " 
      << VolumeNameAndCopy(list->GetDaughter1()) << " " << IsAre(n) << ":\n";

  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetDaught1Points( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }

  out << "    Which in the coordinate system of " 
      << VolumeNameAndCopy(list->GetDaughter2()) << " " << IsAre(n) << ":\n";
  
  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetDaught2Points( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }
  
  out << std::endl;
}


//
// OvershootingDaughter
//
void G4GeomTestStreamLogger::
OvershootingDaughter( const G4GeomTestOvershootList *list )
{
  G4int n = list->NumError();
  if (n <= 0) return;

  out << "GeomTest Error: Overshooting daughter volume\n"
      << "    The volume " << VolumeNameAndCopy(list->GetDaughter())
      << " appears to extend outside the mother volume " 
      << VolumeNameAndCopy(list->GetMother()) << "\n"
      << "    at the following " << (n>1 ? "points" : "point")
      << " in global coordinates:";


  G4int nInterval, nStop;

  if (n <= maxPointsPerError) {
    out << "\n";
    nInterval = 1;
    nStop = n;
  }
  else {
    out << " (list truncated)\n";
    nInterval = n/maxPointsPerError;
    nStop = maxPointsPerError*nInterval;
  }

  G4int i;
  G4ThreeVector s1, s2;
  
  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetGlobalPoints( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }

  out << "    Which in the mother coordinate system " << IsAre(n) << ":\n";

  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetMotherPoints( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }

  out << "    Which in the coordinate system of " 
      << VolumeNameAndCopy(list->GetDaughter()) << " " << IsAre(n) << ":\n";

  PrintSegmentListHeader();
  for(i=0;i<nStop;i+=nInterval) {
    list->GetDaughtPoints( i, s1, s2 );
    PrintSegmentListElement( s1, s2 );
  }
  
  out << std::endl;
}



//
// PrintSegmentListHeader (protected)
//
// Print out a header for a segment list
//
void G4GeomTestStreamLogger::PrintSegmentListHeader()
{
  static const char *header =
     "       length (cm)    ---------- start position (cm) -----------   ----------- end position (cm) ------------\n";
  //      .............|    .............|.............|.............|   .............|.............|.............|
  //  1234              1234                                          123

  out << header;
}


//
// PrintSegmentListElement (protected)
//
// Print out one segment value
//
void G4GeomTestStreamLogger::PrintSegmentListElement( const G4ThreeVector &s1,
                                                      const G4ThreeVector &s2 )
{
  out << "    " << std::setprecision(6) << std::setw(14)
      << (s1-s2).mag()/cm
      << "    " << PrintPos(s1,false) << "   " << PrintPos(s2,false) << "\n";
}


//
// IsAre (protected)
//
// Return a pointer to the string "is" if the argument
// is equal to 1, otherwise return a pointer to "are"
//
const char *G4GeomTestStreamLogger::IsAre( G4int n )
{
  const char *is = "is";
  const char *are = "are";
  
  return n > 1 ? are : is;
}
