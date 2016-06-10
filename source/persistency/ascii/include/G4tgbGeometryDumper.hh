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
// class G4tgbGeometryDumper
//
// Class description:
//
// Class for dumping the whole geometry.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbGeometryDumper_HH
#define G4tgbGeometryDumper_HH

#include "globals.hh"
#include "G4RotationMatrix.hh"

#include <fstream>
#include <map>
#include <vector>

class G4Material;
class G4Element;
class G4Isotope;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVParameterised;
class G4PVReplica;

class G4tgbGeometryDumper
{

  public:  // with description

    static G4tgbGeometryDumper* GetInstance();
    ~G4tgbGeometryDumper();

    void DumpGeometry(const G4String& fname );
    G4VPhysicalVolume* GetTopPhysVol();
    void DumpPhysVol( G4VPhysicalVolume* pv );
    void DumpPVPlacement( G4VPhysicalVolume* pv, const G4String& lvName,
                          G4int copyNo = -999 );
    void DumpPVParameterised( G4PVParameterised* pv );
    void DumpPVReplica( G4PVReplica* pv, const G4String& lvName );
    G4String DumpLogVol( G4LogicalVolume* lv, G4String extraName = "",
                         G4VSolid* solid = 0, G4Material* mate = 0);
    G4String DumpMaterial( G4Material* mat );
    void DumpElement( G4Element* ele);
    void DumpIsotope( G4Isotope* ele);
    G4String DumpSolid( G4VSolid* solid, const G4String& extraName = "" );
    void DumpBooleanVolume( const G4String& solidType, G4VSolid* so );
    void DumpSolidParams(G4VSolid * so);
    std::vector<G4double> GetSolidParams( const G4VSolid * so);
    void DumpPolySections(G4int zPlanes, G4double* z,
                          G4double *rmin, G4double *rmax);
    G4String DumpRotationMatrix( G4RotationMatrix* rotm );

  private:

    G4tgbGeometryDumper();

  private:

    std::vector<G4VPhysicalVolume*> GetPVChildren( G4LogicalVolume* lv );
    G4String GetTGSolidType( const G4String& solidtype );
    G4double MatDeterminant(G4RotationMatrix * ro) ;
    G4double approxTo0( G4double val );
    G4String AddQuotes( const G4String& str );

    G4String GetIsotopeName( G4Isotope* );
    template< class TYP> G4String GetObjectName( TYP* obj,
                         std::map<G4String,TYP*> objectsDumped );
    G4bool CheckIfLogVolExists( const G4String& name, G4LogicalVolume* pt );
    G4bool CheckIfPhysVolExists( const G4String& name, G4VPhysicalVolume* );
    G4String LookForExistingRotation( const G4RotationMatrix* rotm );
    G4String SupressRefl( G4String name );
    G4String SubstituteRefl( G4String name );
    G4bool Same2G4Isotopes( G4Isotope* ele1, G4Isotope* ele2 );
    const G4String& FindSolidName( G4VSolid* solid );

  private:

    static G4ThreadLocal G4tgbGeometryDumper* theInstance;

    std::ofstream* theFile;

    std::map<G4String,G4Material*> theMaterials;
    std::map<G4String,G4Element*> theElements;
    std::map<G4String,G4Isotope*> theIsotopes;
    std::map<G4String,G4VSolid*> theSolids;
    std::map<G4String,G4LogicalVolume*> theLogVols;
    std::map<G4String,G4VPhysicalVolume*> thePhysVols;
    std::map<G4String,G4RotationMatrix*> theRotMats;

    G4int theRotationNumber;
};

#endif
