#ifndef G4tgrFileReader_H
#define G4tgrFileReader_H 1

using namespace std;
#include "globals.hh"
#include <vector>

class G4tgrVolume;

//------------------------------------------------------------------
//
// ClassName:   G4tgrFileReader
//  
// Description: This service provides access to detector description data
//
// Author:      Pedro Arce
// Changes:     17-09-00 Create class
//
//------------------------------------------------------------------

class G4tgrVolumeMgr;

class G4tgrFileReader {

public:

  /// Get the only instance 
  static G4tgrFileReader* GetInstance();  
 
  G4bool Initialize();

  G4tgrFileReader();
  virtual ~G4tgrFileReader(){ };
  void AddTextFile( const G4String& fname ) {
    theTextFiles.push_back( fname );
  } 

private:
  virtual G4bool ProcessLine( const std::vector<G4String>& wl );
  G4tgrVolume* FindVolume( const G4String& volname );

private:
  static G4tgrFileReader* theInstance;  

  vector<G4String> theTextFiles;

  G4tgrVolumeMgr* volmgr;

};


# endif 
