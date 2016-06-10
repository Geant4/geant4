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
// $Id: G4tgrVolumeMgr.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrVolumeMgr
//
// Class description:
//
// Class to manage the detector units. It is a singleton.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrVolumeMgr_h
#define G4tgrVolumeMgr_h

#include "globals.hh"
#include "G4tgrSolid.hh"
#include "G4tgrVolume.hh"
#include "G4tgrPlace.hh"
#include "G4tgrIsotope.hh"
#include "G4tgrElement.hh"
#include "G4tgrMaterial.hh"
#include "G4tgrRotationMatrix.hh"

#include <map>

typedef std::map< G4String, G4tgrSolid* > G4mapssol;
typedef std::map< G4String, G4tgrVolume* > G4mapsvol;
typedef std::multimap< G4String, const G4tgrPlace* > G4mmapspl;

//----------------------------------------------------------------------------  
class G4tgrVolumeMgr 
{ 
  public:  // with description  

    static G4tgrVolumeMgr* GetInstance();  
      // Get the only instance 

    G4tgrSolid* CreateSolid( const std::vector<G4String>& wl, G4bool bVOLUtag );

    void RegisterParentChild( const G4String& parentName,
                              const G4tgrPlace* pl );
      // Add to theG4tgrVolumeTree

    G4tgrSolid* FindSolid( const G4String& name, G4bool exists = false );
      // Find a G4tgrSolid with name 'name'. If it is not found:
      // if exists is true, exit; if exists is false, return 0

    G4tgrVolume* FindVolume( const G4String& volname, G4bool exists = false );
      // Find a G4tgrVolume with name 'volname'. If it is not found:
      // if exists is true, exit; if exists is false, return 0

    std::vector<G4tgrVolume*> FindVolumes( const G4String& volname,
                                                 G4bool exists ); 
      // Find all G4tgrVolume's with name 'volname'. '*' can be used in the 
      // name to mean 'any character' or 'any substring'. If it is not found:
      // if exists is true, exit; if exists is false, return 0

    const G4tgrVolume* GetTopVolume();  
      // Find the top of the volume tree

    std::pair<G4mmapspl::iterator, G4mmapspl::iterator>
    GetChildren( const G4String& name );
      // Find the list of G4tgrPlace children of G4tgrVolume 'name'

    void DumpSummary();
      // Dump summary
    void DumpVolumeTree();
      // Dump to cout the tree of G4tgrVolume's
    void DumpVolumeLeaf( const G4tgrVolume* vol, unsigned int copyNo,
                         unsigned int leafDepth);
      // Dump a G4tgrVolume indicating its copy no
      // and its depth (number of ancestors)

    void RegisterMe( G4tgrSolid* vol);
    void UnRegisterMe( G4tgrSolid* vol );
    void RegisterMe( G4tgrVolume* vol);
    void UnRegisterMe( G4tgrVolume* vol );
    void RegisterMe( G4tgrPlace* pl ) { theG4tgrPlaceList.push_back( pl ); }
    void RegisterMe( G4tgrIsotope* iso ) { theHgIsotList.push_back( iso ); }
    void RegisterMe( G4tgrElement* ele ) { theHgElemList.push_back( ele ); }
    void RegisterMe( G4tgrMaterial* mat ) { theHgMateList.push_back( mat ); }
    void RegisterMe( G4tgrRotationMatrix* rm ) { theHgRotMList.push_back(rm); }

    // Accessors

    const G4mapssol& GetSolidMap() {return theG4tgrSolidMap;}
    const G4mapsvol& GetVolumeMap() {return theG4tgrVolumeMap;}
    const G4mmapspl& GetVolumeTree() {return theG4tgrVolumeTree;}
    std::vector<G4tgrVolume*> GetVolumeList() {return theG4tgrVolumeList;}
    std::vector<G4tgrPlace*> GetDetPlaceList() {return theG4tgrPlaceList;}
    std::vector<G4tgrIsotope*> GetIsotopeList() {return theHgIsotList;}
    std::vector<G4tgrElement*> GetElementList() {return theHgElemList;}
    std::vector<G4tgrMaterial*> GetMaterialList() {return theHgMateList;}
    std::vector<G4tgrRotationMatrix*> GetRotMList() {return theHgRotMList;}

  private:

    G4tgrVolumeMgr();
   ~G4tgrVolumeMgr();

  private:

    G4mapssol theG4tgrSolidMap;
      // Map of G4tgrSolid's: G4String is the G4tgrSolid name,
      // G4tgrSolid* the pointer to it

    G4mapsvol theG4tgrVolumeMap;
      // Map of G4tgrVolume's: G4String is the G4tgrVolume name,
      // G4tgrVolume* the pointer to it

    G4mmapspl theG4tgrVolumeTree;
      // Hierarchy tree of G4tgrVolume's: G4String is the name
      // of the parent G4tgrVolume, G4tgrPlace* the pointers to children

    static G4ThreadLocal G4tgrVolumeMgr* theInstance;

    std::vector<G4tgrVolume*> theG4tgrVolumeList;
    std::vector<G4tgrPlace*> theG4tgrPlaceList;
    std::vector<G4tgrIsotope*> theHgIsotList;
    std::vector<G4tgrElement*> theHgElemList;
    std::vector<G4tgrMaterial*> theHgMateList;
    std::vector<G4tgrRotationMatrix*> theHgRotMList;
};

#endif
