//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// Class description:
//
// Virtual mother class of the division testers for each solid

// History:
// 13.08.03 - P.Arce Initial version
// ********************************************************************

#ifndef ExVDivTester_H
#define ExVDivTester_H 1

#include <vector>
#include <fstream>
#include <iostream>
#include "globals.hh"
#include "geomdefs.hh"

class G4VSolid;
class G4Material;

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VDivisionParameterisation.hh"

enum SolidType{g4box, g4tube, g4tubs, g4cone, g4cons, g4trd, g4para, g4polycone, g4polyhedra };
enum PVType{pvDivision, pvReplica, pvPlacement };
enum PlaceType{pvNormal, pvReflected};
//-enum DivType{divNDiv, divWidth, divNDivWidth };

class ExVDivTester
{ 
  public:  // with description
  
    ExVDivTester( PVType& pvtype, PlaceType& postype,
                  std::vector<G4String>& extraPars );
    virtual ~ExVDivTester();

    void readExtraPars();
    G4double getValueFromExtraPar( G4String& epar );

    void  GenerateRandomPoints();
    virtual void GenerateScanPoints() = 0;
    void GenerateScanPointsAsBox();
    void GenerateScanPointsAsTubs();

    G4VPhysicalVolume* BuildGeometry();
    virtual void BuildParentSolids() = 0;
    virtual void BuildChildrenSolids() = 0;

    void PrintChildrenSolids( std::ostream& out );
    void PrintParentSolid( std::ostream& out );

    inline G4double GetWorldLengthXY() { return theWorldLengthXY; }
    inline G4double GetWorldLengthZ() { return theWorldLengthZ; }
    inline G4double GetWorldGap() { return theWorldGap; }

  private:
    void BuildParentVolumes( G4LogicalVolume* worldLog );
    void BuildChildrenVolumes();
    void SetOffsets();

  public:
    static G4bool bDivCylindrical;

  protected:

    PVType thePVType;
    PlaceType thePosType;
    DivisionType theDivType;
    std::vector<G4String> theExtraPars;
    std::vector<G4double> theWidths;
    std::vector<G4double> theOffsets;
    G4double theOffsetFactor;
    G4double theWidthFactor;
    G4int theNDiv;
    G4double theStartPhi;
    G4double theDeltaPhi;

    std::vector<G4VSolid*> theParentSolids;
    std::vector<G4LogicalVolume*> theParentLogs;
    std::vector<G4VPhysicalVolume*> theParentPhyss;
    std::vector<G4VSolid*> theChildSolids;
    std::vector<G4LogicalVolume*> theChildLogs;
    std::vector<EAxis> theAxis;

    G4double theWorldLengthXY;
    G4double theWorldLengthZ;
    G4double theWorldGap;

    G4Material* theMate;

    G4int numberOfPoints;

    G4int verbose;
};

#endif

