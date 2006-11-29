// Include files
#include "G4tgrFileReader.hh"
#include "G4tgrParameterMgr.hh"
#include "G4tgrFileIn.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrVolume.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrVolumeDivision.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"

G4tgrFileReader* G4tgrFileReader::theInstance = 0;


//---------------------------------------------------------------
G4tgrFileReader::G4tgrFileReader()
{
  new G4tgrMessenger;
}


//---------------------------------------------------------------
G4tgrFileReader* G4tgrFileReader::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgrFileReader;
  }
  return theInstance;
}

//---------------------------------------------------------------
bool G4tgrFileReader::Initialize() 
{

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgrFileReader::Initialize " << G4endl;
#endif

  volmgr = G4tgrVolumeMgr::GetInstance();
  vector< G4String > wl,wlnew;
    
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " Number of geometry data files = " << theTextFiles.size() << G4endl;
#endif

  for( uint ii = 0; ii < theTextFiles.size(); ii++ ) {
    
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 ) {
      G4cout << " Data File is " << theTextFiles[ii] << G4endl;
      //      G4cout << "We should open the file here" << G4endl;
    }
#endif
    
    G4tgrFileIn fin = G4tgrFileIn::GetInstance( theTextFiles[ii] );
    
    int nlines = 0;
    for(;;) {
      nlines++;
      if(! fin.GetWordsInLine( wlnew ) ) break; 
      // Check if it is continuation line or first line
      if( wlnew[0].c_str()[0] != ':' ) {
        wl.insert( wl.end(), wlnew.begin(), wlnew.end() );
	if( G4tgrMessenger::GetVerboseLevel() >= 3 ) G4tgrUtils::DumpVS(wl, "!!!! adding line");
	if( G4tgrMessenger::GetVerboseLevel() >= 3 ) G4tgrUtils::DumpVS(wlnew, "!!!! line new");
	continue;
      }else {
	//----- Process previous tag
	if( G4tgrMessenger::GetVerboseLevel() >= 2 ) G4tgrUtils::DumpVS(wl, "!!!! line read");
	if( nlines != 1) { // first line has no previous tag
	  if( ! ProcessLine( wl ) ) {
	    fin.DumpException( "tag not found" );
	  }
	}
	wl = wlnew;
      }
    
    }
    
    if( wl.size() != 0 ) {
      if( ! ProcessLine( wl ) ) {
	fin.DumpException( "tag not found" );
      }
    }

  }
  
  return 1;
}


//---------------------------------------------------------------
bool G4tgrFileReader::ProcessLine( const std::vector<G4String>& wl )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4tgrUtils::DumpVS(wl, "!!!! processing line");
#endif  

  G4String wl0 = wl[0];
  for( size_t ii = 0; ii < wl0.length(); ii++ ){
    wl0[ii] = toupper( wl0[ii] );
  }

  //------------------------------- parameter
  if( wl0 == ":P" ) {
    G4tgrParameterMgr::GetInstance()->AddParameter( wl );
    
  //------------------------------- isotope
  }else if( wl0 == ":ISOT" ) {
    G4tgrIsotope* isot = G4tgrMaterialFactory::GetInstance()->AddIsotope( wl );
    volmgr->RegisterMe( isot );

  //------------------------------- element
  }else if( wl0 == ":ELEM" ) {
    G4tgrElementSimple* elem = G4tgrMaterialFactory::GetInstance()->AddElementSimple( wl );
    volmgr->RegisterMe( elem );

  //------------------------------- element from isotopes
  }else if( wl0 == ":ELEM_FROM_ISOT" ) {
    //:ELEM_FROM_ISOT NAME SYMBOL N_ISOT (ISOT_NAME ISOT_ABUNDANCE)
    G4tgrElementFromIsotopes* elem = G4tgrMaterialFactory::GetInstance()->AddElementFromIsotopes( wl );
    volmgr->RegisterMe( elem );
    
  //------------------------------- material
  } else if( wl0 == ":MATE" ) {
    G4tgrMaterialSimple* mate = G4tgrMaterialFactory::GetInstance()->AddMaterialSimple( wl );
    volmgr->RegisterMe( mate );
    
    //------------------------------- material mixtures
    //------------------------------- material mixture by weight
  } else if( wl0 == ":MIXT" 
	     || wl0 == ":MIXT_BY_WEIGHT" ) {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
     G4tgrUtils::DumpVS( wl, " material mixture by weight read " );
#endif
    G4tgrMaterialMixture* mate = G4tgrMaterialFactory::GetInstance()->AddMaterialMixture( wl, "MaterialMixtureByWeight" );
    volmgr->RegisterMe( mate );

    //------------------------------- material mixture by number of atoms
  } else if( wl0 == ":MIXT_BY_NATOMS" ) {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4tgrUtils::DumpVS( wl, " material mixture by number of atoms read " );
#endif
    G4tgrMaterialMixture* mate = G4tgrMaterialFactory::GetInstance()->AddMaterialMixture( wl, "MaterialMixtureByNoAtoms" );
    volmgr->RegisterMe( mate );

    //------------------------------- material mixture by volume
  } else if( wl0 == ":MIXT_BY_VOLUME" ) {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4tgrUtils::DumpVS( wl, " material mixture by volume read " );
#endif
    G4tgrMaterialMixture* mate = G4tgrMaterialFactory::GetInstance()->AddMaterialMixture( wl, "MaterialMixtureByVolume" );
    volmgr->RegisterMe( mate );

    //------------------------------- solid
  } else if( wl0 == ":SOLID" ) {

    G4tgrSolid* sol = volmgr->CreateSolid( wl, 0 ); // called from here or from G4tgrVolume::G4tgrVolume

    //------------------------------- volume
  } else if( wl0 == ":VOLU" ) {

    G4tgrVolume* vol = new G4tgrVolume( wl );

    volmgr->RegisterMe( vol );
    
    //--------------------------------- single placement
  } else if( wl0 == ":PLACE" ) {
    G4tgrVolume* vol = FindVolume( G4tgrUtils::SubQuotes( wl[1] ) );
    G4tgrPlace* vpl = vol->AddPlace( wl );
    //t ?????    
    volmgr->RegisterMe( vpl );
    
    //--------------------------------- parameterisation
  } else if( wl0 == ":PLACE_PARAM" ) {
    G4tgrVolume* vol = FindVolume( G4tgrUtils::SubQuotes( wl[1] ) );
    G4tgrPlaceParameterisation* vpl = vol->AddPlaceParam( wl );
    volmgr->RegisterMe( vpl );
	
    //--------------------------------- division
  } else if( wl0 == ":DIV_NDIV" 
	     || wl0 == ":DIV_WIDTH" 
	     || wl0 == ":DIV_NDIV_WIDTH" ) {
    //---------- Create G4tgrVolumeDivision and fill the volume params
    G4tgrVolumeDivision* vol = new G4tgrVolumeDivision( wl );
    volmgr->RegisterMe( vol );

    //--------------------------------- replica
  } else if( wl0 == ":REPL" ) {
    G4tgrVolume* vol = FindVolume( G4tgrUtils::SubQuotes( wl[1] ) );
    G4tgrPlaceDivRep* vpl = vol->AddPlaceReplica( wl );
    //    G4tgrPlaceDivRep* vpl = new G4tgrPlaceDivRep( wl );
    volmgr->RegisterMe( vpl );
    
    //---------------------------------  rotation matrix
  } else if( wl0 == ":ROTM" ) {
    //---------- When second word is ':NEXT/:MNXT' it is used for defining a rotation matrix that will be used for the next placement/s
    G4tgrRotationMatrix* rm = G4tgrRotationMatrixFactory::GetInstance()->AddRotMatrix( wl );
    volmgr->RegisterMe( rm );
    
    //------------------------------- visualisation
  } else if( wl0 == ":VIS" ) {
    G4tgrVolume* vol = volmgr->FindVolume( G4tgrUtils::SubQuotes( wl[1] ), 1 );
    vol->AddVisibility( wl );

    //--------------------------------- colour
  } else if( wl0 == ":COLOUR"|| wl0 == ":COLOR" ) {
    G4tgrVolume* vol = volmgr->FindVolume( G4tgrUtils::SubQuotes( wl[1] ), 1 );
    vol->AddRGBColour( wl );
    
    //--------------------------------- ERROR
  } else {
    G4tgrUtils::DumpVS( wl, "!!! G4TgFileReader TAG NOT SUPPORTED ");
    return 0;
  } 

  return 1;
  
}

//---------------------------------------------------------------
G4tgrVolume* G4tgrFileReader::FindVolume( const G4String& volname )
{
  G4tgrVolume* vol;

  G4tgrVolume* volt = volmgr->FindVolume( volname, 1);

  if( volt->GetType() == "VOLDivision" ) {
    G4Exception(" !!!EXITING: G4tgrFileReader:initialize. Using 'PLACE' for a volume created by a division ");
  } else { 
    vol = (G4tgrVolume*)volt;
  }

  return vol;
}



