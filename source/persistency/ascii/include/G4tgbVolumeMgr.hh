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
// $Id: G4tgbVolumeMgr.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgbVolume
//
// Class description:
//
// Class to manage volumes: G4VSolids, G4LogicalVolumes, G4VPhysicalVolumes.
// It is a singleton, accesed always through calls to GetInstance().

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbVolumeMgr_h
#define G4tgbVolumeMgr_h

#include "globals.hh"

#include <string>
#include <vector>
#include <map>

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

class G4tgbVolume;
class G4tgrVolume;
class G4tgbDetectorBuilder;

typedef std::map< G4String, G4tgbVolume* > G4mssvol;
typedef std::multimap< G4String, G4VSolid* > G4mmssol;
typedef std::multimap< G4String, G4LogicalVolume* > G4mmslv;
typedef std::multimap< G4String, G4VPhysicalVolume* > G4mmspv;
typedef std::map< G4LogicalVolume*, G4LogicalVolume* > G4mlvlv;
typedef std::map< G4VPhysicalVolume*, G4VPhysicalVolume* > G4mpvpv;

//----------------------------------------------------------------------------  
class G4tgbVolumeMgr 
{ 
  public:  // with description

    G4tgbVolumeMgr();
   ~G4tgbVolumeMgr();

    static G4tgbVolumeMgr* GetInstance();  
      // Get the only instance 

    void AddTextFile( const G4String& fname );
    G4VPhysicalVolume* ReadAndConstructDetector();

    void CopyVolumes();
      // Build a G4tgbVolume per each G4tgbVolume

    G4tgbVolume* FindVolume( const G4String& volname);
      // Find a G4tgbVolume by name

    void RegisterMe( const G4tgbVolume* vol );
      // Register a G4tgbVolume
    void RegisterMe( const G4VSolid* solid );
      // Register a G4VSolid
    void RegisterMe( const G4LogicalVolume* lv );
      // Register a G4LogicalVolume
    void RegisterMe( const G4VPhysicalVolume* pv );
      // Register a G4VPhysicalVolume
    void RegisterChildParentLVs( const G4LogicalVolume* logvol,
                                 const G4LogicalVolume* parentLV );
      // Register a child and its parent LV

    G4VSolid* FindG4Solid( const G4String& name );
      // Find if solid already exists, comparing the name and all parameters 
      // (could be checked before creating it, but it would be quite
      // complicated, because it would have to compare the parameters, and
      // they depend on the type of solid)

    G4LogicalVolume* FindG4LogVol( const G4String& theName,
                                   const G4bool bExists = 0 );
      // Find a G4LogicalVolume if it already exists

    G4VPhysicalVolume* FindG4PhysVol( const G4String& theName,
                                   const G4bool bExists = 0 );
      // Find a G4VPhysicalVolume if it already exists

    G4VPhysicalVolume* GetTopPhysVol();
      // Get the top PV in the hierarchy tree: calls topLV, because
      // physicalvolumes are not placed until geometry is initialized

    G4LogicalVolume* GetTopLogVol();
      // Get the top LV in the hierarchy tree

    void BuildPhysVolTree();

    // Dumping methods

    void DumpSummary();
    void DumpG4LogVolTree();
    void DumpG4LogVolLeaf(const G4LogicalVolume* lv, unsigned int leafDepth);
    void DumpG4PhysVolTree();
    void DumpG4PhysVolLeaf(const G4VPhysicalVolume* pv, unsigned int leafDepth);
    void DumpG4SolidList();

  public:  // without description

    const std::multimap< G4String, G4VSolid* >& GetSolids() const
      { return theSolids; }
    void SetDetectorBuilder( G4tgbDetectorBuilder* db )
      { theDetectorBuilder = db; }
    G4tgbDetectorBuilder* GetDetectorBuilder() const
      { return theDetectorBuilder; }

  private:

    static G4ThreadLocal G4tgbVolumeMgr* theInstance;

    G4mssvol theVolumeList;
      // Map of G4tgbVolume's: G4String is the G4tgbVolume name,
      // G4tgbVolume* the pointer to it.

    G4mmssol theSolids;
      // Solids container

    G4mmslv theLVs;
      // Logical volume container
    G4mmspv thePVs;
      // Physical volume container

    G4mlvlv theLVTree;
      // Logical volume tree for navigation (from parent to children):
      // first is parent, then child
    G4mlvlv theLVInvTree;
      // Logical volume tree for inverse navigation (from children to parent):
      // first is child, then parent

    G4mpvpv thePVTree;
      // Physical volume tree for navigation (from parent to children):
      // first is parent, then child
    G4mpvpv thePVInvTree;
      // Physical volume tree for inverse navigation (from children to parents):
      // first is child, then parent

    G4tgbDetectorBuilder* theDetectorBuilder;
};

#endif
