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
#include "DicomROIContour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomROIContour::DicomROIContour()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomROIContour::SetData( std::vector<G4ThreeVector> data )
{
  thePoints = data;
  if( theGeomType == "CLOSED_PLANAR" ){
    //add a last point to close the circle
    thePoints.push_back( thePoints[0] );
    for( size_t ii = 1; ii < thePoints.size(); ii++ ){
      theDirections.push_back( thePoints[ii] - thePoints[ii-1] );
    }
    theDirections.push_back( thePoints[0] - thePoints[thePoints.size()-1] ); // (0,0)
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DicomROIContour::GetZ()
{
  if( thePoints.size() == 0 ) {
    G4Exception("DicomROIContour::GetZ",
                "DRC001",
                FatalException,
                "No Point in contour");
    Print(G4cout);
    return DBL_MAX;
  }
  return thePoints[0].z();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomROIContour::AddPoints( std::vector<G4ThreeVector> points )
{
  points.push_back( points[0] );
  size_t npold = thePoints.size();
  for( size_t ii = 0; ii < points.size(); ii++ ) {
    thePoints.push_back( points[ii] );
  }

  for( size_t ii = npold+1; ii < thePoints.size(); ii++ ){
    theDirections.push_back( thePoints[ii] - thePoints[ii-1] );
  }
  theDirections.push_back( thePoints[npold] - thePoints[thePoints.size()-1] );
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomROIContour::Print( std::ostream& out )
{

  out << this << "@@@ NUMBER OF ContourImageSequences " << theImageIUIDs.size() << G4endl;
  for( size_t ii = 0; ii < theImageIUIDs.size(); ii++ ) {
    out << "ContourImageSequence I= " << theImageIUIDs[ii] << G4endl;
  }
  out << "@@@ GeomType " << theGeomType << " NPOINTS " << thePoints.size() << G4endl;
  for( size_t ii = 0; ii < thePoints.size(); ii++ ) {
    out << "(" << thePoints[ii].x()<< "," << thePoints[ii].y()<< "," << thePoints[ii].z() 
        << ")" << G4endl;
  }

}
