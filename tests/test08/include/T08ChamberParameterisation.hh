// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08ChamberParameterisation.hh,v 1.1 1999-01-08 16:35:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal width, & their lengths are a linear equation.
//    They are spaced an equal distance apart, starting from given location.
//
#ifndef T08ChamberParameterisation_H
#define T08ChamberParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
class G4VPhysicalVolume;
class G4Box;

class T08ChamberParameterisation : public G4VPVParameterisation
{ 
  public:
    T08ChamberParameterisation(G4int    NoChambers, 
                                 G4double startZ, 
                                 G4double spacing,
                                 G4double widthChamber, 
                                 G4double lengthInitial,
                                 G4double lengthFinal );
    ~T08ChamberParameterisation();
    void ComputeTransformation
    (const G4int copyNo,G4VPhysicalVolume *physVol) const;
    void ComputeDimensions
    (G4Box & trackerLayer, const G4int copyNo,
      const G4VPhysicalVolume * physVol) const;

    // Functions to get the parameters would be nice.
  // eg  G4int GetNoChambers();

  private:

    G4double fNoChambers;   //  
    G4double fStartZ;
    G4double fHalfWidth;    //  The half-width of each tracker chamber
    G4double fSpacing;      //  The distance between the chambers' center
    G4double fHalfLengthFirst;  //  The first half-length 
    G4double fHalfLengthIncr;   //  The Increment for the half-length 
};

#endif


