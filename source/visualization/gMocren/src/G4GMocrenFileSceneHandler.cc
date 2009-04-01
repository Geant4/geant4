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
// $Id: G4GMocrenFileSceneHandler.cc,v 1.1 2009-04-01 13:16:11 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Akinori Kimura    March 31, 2009
//
// GMocrenFile scene.


//----- header files
#include <fstream>
#include <cstdlib>
#include <cstring>
#include "globals.hh"
#include "G4FRConst.hh"
#include "G4GMocrenFile.hh"
#include "G4GMocrenFileSceneHandler.hh"
#include "G4GMocrenFileViewer.hh"
#include "G4Point3D.hh"
#include "G4VisAttributes.hh"
#include "G4Scene.hh"
#include "G4Transform3D.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyline.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include "G4VPVParameterisation.hh"
#include "G4VVolumeMaterialScanner.hh"
#include "G4VisTrajContext.hh"
#include "G4VisManager.hh"
#include "G4VTrajectoryModel.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "G4HitsModel.hh"
#include "G4GMocrenMessenger.hh"
#include "G4ScoringHitsModel.hh"
#include "G4GMocrenIO.hh"
#include "G4VNestedParameterisation.hh"
#include "G4GMocrenTouchable.hh"

//----- constants
const char  FR_ENV_CULL_INVISIBLE_OBJECTS [] = "G4DAWN_CULL_INVISIBLE_OBJECTS";
const char  GDD_FILE_HEADER      [] = "g4_";
const char  DEFAULT_GDD_FILE_NAME[] = "g4_00.gdd";

const int   FR_MAX_FILE_NUM = 100 ;

const bool GFDEBUG_TRK = false;//true;
const bool GFDEBUG_HIT = false;//true;
//const bool GFDEBUG = false;
const G4int GFDEBUG_DET = 0; // 0: false 

///////////////////////////
// Driver-dependent part //
///////////////////////////


//----- G4GMocrenFileSceneHandler, constructor
G4GMocrenFileSceneHandler::G4GMocrenFileSceneHandler(G4GMocrenFile& system,
						     G4GMocrenMessenger & messenger,
						     const G4String& name)
  : G4VSceneHandler(system, fSceneIdCount++, name) ,
    fSystem(system),
    fMessenger(messenger),
    fgMocrenIO(new G4GMocrenIO()),
    fbSetModalityVoxelSize(false),
    fbModelingTrajectory(false),
    fGddDest(),
    FRflag_in_modeling(false),
    flag_saving_g4_gdd(false),
    COMMAND_BUF_SIZE(G4FRofstream::SEND_BUFMAX),
    fPrec(9), fPrec2(16) {

  // g4.gdd filename and its directory
  if ( getenv( "G4GMocrenFile_DEST_DIR" ) == NULL ) {
    std::strcpy( fGddDestDir , "" )                      ;  // output dir
    std::strcpy( fGddFileName, DEFAULT_GDD_FILE_NAME );  // filename
  } else {
    std::strcpy( fGddDestDir , getenv( "G4GMocrenFile_DEST_DIR" ) ); // output dir
    std::strcpy( fGddFileName, DEFAULT_GDD_FILE_NAME        ); // filename 
  }
		
  // maximum number of g4.gdd files in the dest directory
  fMaxFileNum = FR_MAX_FILE_NUM ; // initialization
  if ( getenv( "G4GMocrenFile_MAX_FILE_NUM" ) != NULL ) {	
		
    std::sscanf( getenv("G4GMocrenFile_MAX_FILE_NUM"), "%d", &fMaxFileNum ) ;

  } else {
    fMaxFileNum = FR_MAX_FILE_NUM ;
  }
  if( fMaxFileNum < 1 ) { fMaxFileNum = 1 ; }

  /*

  //----- precision control
  if( getenv( "G4GMocrenFile_PRECISION" ) != NULL ) {
    sscanf( getenv("G4GMocrenFile_PRECISION"), "%d", &fPrec ) ;
  } else {
    fPrec = 9 ;
  }
  fPrec2 = fPrec + 7 ;
  */
  for(int i = 0; i < 3; i++) {
    fModalitySize[i] = 0;
  }
} 


//----- G4GMocrenFileSceneHandler, destructor
G4GMocrenFileSceneHandler::~G4GMocrenFileSceneHandler () 
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** ~G4GMocrenFileSceneHandler" << G4endl;
#endif 
  if (fGddDest.IsOpen()) 
    {
      //----- End of modeling
      // !EndModeling, !DrawAll, !CloseDevice,
      // close g4.gdd
      FREndModeling();
    }
  ClearStore (); // clear current scene
}

//-----
void	G4GMocrenFileSceneHandler::SetGddFileName() 
{
  // g4_00.gdd, g4_01.gdd, ..., g4_MAX_FILE_INDEX.gdd
  const int MAX_FILE_INDEX = fMaxFileNum - 1 ;

  // dest directory (null if no environmental variables is set)
  std::strcpy ( fGddFileName, fGddDestDir) ; 

  // create full path name (default)
  std::strcat ( fGddFileName, DEFAULT_GDD_FILE_NAME );

  // Automatic updation of file names
  for( int i = 0 ; i < fMaxFileNum ; i++) { 

    // Message in the final execution
    if( i == MAX_FILE_INDEX ) 
      {
	G4cerr << "==========================================="   << G4endl; 
	G4cerr << "WARNING MESSAGE from GMocrenFile driver:   "   << G4endl;
	G4cerr << "  This file name is the final one in the   "   << G4endl;
	G4cerr << "  automatic updation of the output file name." << G4endl; 
	G4cerr << "  You may overwrite existing files, i.e.   "   << G4endl; 
	G4cerr << "  g4_XX.gdd."   << G4endl;
	G4cerr << "==========================================="   << G4endl; 
      }

    // re-determine file name as G4GMocrenFile_DEST_DIR/g4_XX.gdd 
    if( i >=  0 && i <= 9 ) { 
      std::sprintf( fGddFileName, "%s%s%s%d.gdd" , fGddDestDir,  GDD_FILE_HEADER, "0", i );
    } else {
      std::sprintf( fGddFileName, "%s%s%d.gdd" , fGddDestDir,  GDD_FILE_HEADER, i );
    }

    // check validity of the file name
    std::ifstream  fin ; 
    fin.open(fGddFileName) ;
    if(!fin) { 
      // new file	
      fin.close();  
      break; 
    } else { 
      // already exists (try next) 
      fin.close(); 
    } 

  } // for 

  G4cerr << "===========================================    " << G4endl; 
  G4cerr << "Output file: " <<    fGddFileName             << G4endl; 
  G4cerr << "Destination directory (current dir if NULL): "       << fGddDestDir    << G4endl; 
  G4cerr << "Maximal number of files in the destination directory: " << fMaxFileNum << G4endl; 
  G4cerr << "Note:                                                " << G4endl; 
  G4cerr << "  * The maximal number is customizable as:           " << G4endl;
  G4cerr << "       % setenv  G4GMocrenFile_MAX_FILE_NUM  number " << G4endl;        
  G4cerr << "  * The destination directory is customizable as:" << G4endl;
  G4cerr << "       % setenv  G4GMocrenFile_DEST_DIR  dir_name/  " << G4endl;        
  G4cerr << "     ** Do not forget \"/\" at the end of the    " << G4endl;              
  G4cerr << "        dir_name, e.g. \"./tmp/\".  " << G4endl;              
  G4cerr << "===========================================      " << G4endl; 

} // G4GMocrenFileSceneHandler::SetGddFileName()


//-----
void	G4GMocrenFileSceneHandler::BeginSavingGdd( void ) 
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** BeginSavingGdd (called)\n";
#endif

  if( !IsSavingGdd() ) {

#if defined DEBUG_FR_SCENE
    G4cerr << "*****                   (started) " ;
    G4cerr << "(open g4.gdd, ##)"  << G4endl;
#endif
    SetGddFileName() ; // result set to fGddFileName
    fGddDest.Open(fGddFileName)   ;

    //SendStr( FR_G4_GDD_HEADER   )    ; 
    flag_saving_g4_gdd = true        ; 

    G4String fname = "modality-map.dat";
    std::ifstream ifile(fname);
    if(ifile) {
      short minmax[2];
      ifile >> minmax[0] >> minmax[1];
      fgMocrenIO->setModalityImageMinMax(minmax);
      std::vector<float> map;
      float dens;
      for(int i = minmax[0]; i <= minmax[1]; i++) {
	ifile >> dens;
	map.push_back(dens);
      }
      fgMocrenIO->setModalityImageDensityMap(map);
      
    } else {
      G4cerr << "cann't open the file : " << fname << G4endl;
    }

    // mesh size
    //fMessenger.getNoVoxels(fModalitySize[0], fModalitySize[1], fModalitySize[2]);
    //fgMocrenIO->setModalityImageSize(fModalitySize);
    
    // initializations
    fgMocrenIO->clearTracks();
    fgMocrenIO->clearDetector();
    std::vector<Detector>::iterator itr = fDetectors.begin();
    for(; itr != fDetectors.end(); itr++) {
      itr->clear();
    }
    fDetectors.clear();
    
    fNestedHitsList.clear();
    fNestedVolumeNames.clear();
      
  }
}

void	G4GMocrenFileSceneHandler::EndSavingGdd  ( void ) 
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** EndSavingGdd (called)\n";
#endif

  if(  IsSavingGdd() )
    {
#if defined DEBUG_FR_SCENE
      G4cerr << "*****                 (started) (close g4.gdd)" << G4endl;
#endif
      fGddDest.Close()               ;
      flag_saving_g4_gdd = false ; 

      std::map<Index3D, float>::iterator itr = fNestedModality.begin();
      G4int xmax=0, ymax=0, zmax=0;
      for(; itr != fNestedModality.end(); itr++) {
	if(itr->first.x > xmax) xmax = itr->first.x;
	if(itr->first.y > ymax) ymax = itr->first.y;
	if(itr->first.z > zmax) zmax = itr->first.z;
      }
      // mesh size
      fModalitySize[0] = xmax+1;
      fModalitySize[1] = ymax+1;
      fModalitySize[2] = zmax+1;
      fgMocrenIO->setModalityImageSize(fModalitySize);
      G4cout << " modality size : "
	     << fModalitySize[0] << " x "
	     << fModalitySize[1] << " x "
	     << fModalitySize[2] << G4endl;

      G4int nxy = fModalitySize[0]*fModalitySize[1];
      //std::map<G4int, float>::iterator itr;
      for(int z = 0; z < fModalitySize[2]; z++) {
	short * modality = new short[nxy];
	for(int y = 0; y < fModalitySize[1]; y++) {
	  for(int x = 0; x < fModalitySize[0]; x++) {
	  //for(int x = fModalitySize[0]-1; x >= 0 ; x--) {
	    //G4int ixy = x + (fModalitySize[1]-y-1)*fModalitySize[0];

	    G4int ixy = x + y*fModalitySize[0];
	    Index3D idx(x,y,z);
	    itr = fNestedModality.find(idx);
	    if(itr != fNestedModality.end()) {

	      modality[ixy] = fgMocrenIO->convertDensityToHU(itr->second);
	    } else {
	      G4cout << "ABC : " << x << ", " <<  y << ", " << z << G4endl;
	      modality[ixy] = -1024;
	    }

	  }
	}
	fgMocrenIO->setModalityImage(modality);
      }

      //-- dose
      size_t nhits = fNestedHitsList.size();
      G4cout << " # hits : " << nhits << G4endl;

      std::map<Index3D, G4double>::iterator hitsItr;
      std::map<G4String, std::map<Index3D, G4double> >::iterator hitsListItr = fNestedHitsList.begin();

      for(int n = 0; hitsListItr != fNestedHitsList.end(); hitsListItr++, n++) {

	fgMocrenIO->newDoseDist();
	fgMocrenIO->setDoseDistName(hitsListItr->first, n);
	fgMocrenIO->setDoseDistSize(fModalitySize, n);

	G4double minmax[2] = {DBL_MAX, -DBL_MAX};
	for(int z = 0 ; z < fModalitySize[2]; z++) {
	  G4double * values = new G4double[nxy];
	  for(int y = 0; y < fModalitySize[1]; y++) {
	    for(int x = 0; x < fModalitySize[0]; x++) {

	      G4int ixy = x + y*fModalitySize[0];
	      Index3D idx(x,y,z);
	      hitsItr = hitsListItr->second.find(idx);
	      if(hitsItr != hitsListItr->second.end()) {

		values[ixy] = hitsItr->second;
	      } else {
		values[ixy] = 0.;
	      }
	      if(values[ixy] < minmax[0]) minmax[0] = values[ixy];
	      if(values[ixy] > minmax[1]) minmax[1] = values[ixy];
	    }
	  }
	  fgMocrenIO->setDoseDist(values, n);
	}
	fgMocrenIO->setDoseDistMinMax(minmax, n);
	G4double lower = 0.;
	if(minmax[0] < 0)  lower = minmax[0];
	G4double scale = (minmax[1]-lower)/25000.;
	fgMocrenIO->setDoseDistScale(scale, n);
	G4String sunit("unit?"); //temporarily
	fgMocrenIO->setDoseDistUnit(sunit, n);
      }
      

      //-- draw axes
      if(false) {//true) {
      G4ThreeVector trans;
      G4RotationMatrix rot;
      trans = fVolumeTrans3D.getTranslation();
      rot = fVolumeTrans3D.getRotation().inverse();
      // x
      std::vector<float *> tracks;
      unsigned char colors[3];
      float * trk = new float[6];
      tracks.push_back(trk);

      G4ThreeVector orig(0.,0.,0), xa(1000.,0.,0.), ya(0.,1000.,0.), za(0.,0.,1000.);
      orig -= trans;
      orig.transform(rot);
      xa -= trans;
      xa.transform(rot);
      ya -= trans;
      ya.transform(rot);
      za -= trans;
      za.transform(rot);
      for(int i = 0; i < 3; i++) trk[i] = orig[i];
      for(int i = 0; i < 3; i++) trk[i+3] = xa[i];
      colors[0] = 255; colors[1] = 0; colors[2] = 0;
      fgMocrenIO->addTrack(tracks, colors);
      // y
      for(int i = 0; i < 3; i++) trk[i+3] = ya[i];
      colors[0] = 0; colors[1] = 255; colors[2] = 0;
      fgMocrenIO->addTrack(tracks, colors);
      // z
      for(int i = 0; i < 3; i++) trk[i+3] = za[i];
      colors[0] = 0; colors[1] = 0; colors[2] = 255;
      fgMocrenIO->addTrack(tracks, colors);
      }

      //-- detector
      ExtractDetector();


      if(GFDEBUG_DET) G4cout << ">>>>>>>>>>>>>>>>>>>>>>   (";
      std::vector<float> transformObjects;
      for(int i = 0; i < 3; i++) {
	// need to check!!
	transformObjects.push_back((fVolumeSize[i]/2. - fVoxelDimension[i]/2.));
	if(GFDEBUG_DET) G4cout << transformObjects[i] << ", ";
      }
      if(GFDEBUG_DET) G4cout << ")" << G4endl;


      fgMocrenIO->translateTracks(transformObjects);
      fgMocrenIO->translateDetector(transformObjects);

      // store
      fgMocrenIO->storeData(fGddFileName);
    } 
}


//----- 
void G4GMocrenFileSceneHandler::FRBeginModeling( void )
{
  G4VSceneHandler::BeginModeling();

  if( !FRIsInModeling() )  	
    {
#if defined DEBUG_FR_SCENE
      G4cerr << "***** G4GMocrenFileSceneHandler::FRBeginModeling (called & started)" << G4endl;
#endif


      //----- Send saving command and heading comment
      BeginSavingGdd();

      FRflag_in_modeling = true ;


      G4Scene * scene = GetScene();
      //G4VModel * model = new G4ScoringHitsModel();
      //scene->AddEndOfEventModel(model);
      scene->AddEndOfEventModel(new G4ScoringHitsModel());
      scene->AddEndOfEventModel(new G4HitsModel());
      if(GFDEBUG_HIT) {
	std::vector<G4VModel*> vmodel = scene->GetEndOfEventModelList();
	std::vector<G4VModel*>::iterator itr = vmodel.begin();
	for(; itr != vmodel.end(); itr++) {
	  G4cout << " IIIIII model name: " << (*itr)->GetGlobalTag() << G4endl;
	}
      }

    } // if


} 


//========== AddPrimitive() functions ==========//

//----- Add polyline 
void G4GMocrenFileSceneHandler::AddPrimitive (const G4Polyline& polyline) 
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive\n";
#endif 

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();


  // draw trajectories
  if(fbModelingTrajectory) {
    
    G4TrajectoriesModel * pTrModel = dynamic_cast<G4TrajectoriesModel*>(fpModel);
    if (!pTrModel) { 
      G4Exception 
	("G4VSceneHandler::AddCompound(const G4Polyline&): Not a G4TrajectoriesModel.");
    }

    G4ThreeVector trans;
    G4RotationMatrix rot;
    trans = fVolumeTrans3D.getTranslation();
    rot = fVolumeTrans3D.getRotation().inverse();

    if(GFDEBUG_TRK) G4cout << "   trajectory points : " << G4endl;
    std::vector<float *> trajectory;
    if(polyline.size() < 2) return;
    G4Polyline::const_iterator preitr = polyline.begin();
    G4Polyline::const_iterator postitr = preitr; postitr++;
    for(; postitr != polyline.end(); preitr++, postitr++) {
      G4ThreeVector prePts(preitr->x(), preitr->y(), preitr->z());
      prePts -= trans;
      prePts.transform(rot);
      G4ThreeVector postPts(postitr->x(), postitr->y(), postitr->z());
      postPts -= trans;
      postPts.transform(rot);
      float * stepPts = new float[6];
      stepPts[0] = prePts.x();
      stepPts[1] = prePts.y();
      stepPts[2] = prePts.z();
      stepPts[3] = postPts.x();
      stepPts[4] = postPts.y();
      stepPts[5] = postPts.z();
      trajectory.push_back(stepPts);

      if(GFDEBUG_TRK) {
	G4cout << "                       ("
	       << stepPts[0] << ", "
	       << stepPts[1] << ", "
	       << stepPts[2] << ") - (" 
	       << stepPts[3] << ", "
	       << stepPts[4] << ", "
	       << stepPts[5] << ")" << G4endl;
      }
    }

    const G4VisAttributes * att = polyline.GetVisAttributes();
    G4Color color = att->GetColor();
    unsigned char trkcolor[3];
    trkcolor[0] = (unsigned char)(color.GetRed()*255);
    trkcolor[1] = (unsigned char)(color.GetGreen()*255);
    trkcolor[2] = (unsigned char)(color.GetBlue()*255);
    if(GFDEBUG_TRK) {
      G4cout << "              color  : ["
	     << color.GetRed() << ", "
	     << color.GetGreen() << ", "
	     << color.GetBlue() << "]" << G4endl;
    }

    fgMocrenIO->addTrack(trajectory, trkcolor);

  }

} // G4GMocrenFileSceneHandler::AddPrimitive (polyline)


//----- Add nurbes
void G4GMocrenFileSceneHandler::AddPrimitive (const G4NURBS&)
{
  //----- 
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive( G4NURBS )\n";
#endif

  //----- Initialize DAWN IF NECESSARY
  FRBeginModeling();
	
  ///////////////////////////////////////////////
  // DAWN does not support NUBS visualizaition //
  ///////////////////////////////////////////////
}



//----- Add text
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Text& text )
{
  // to avoid a warning in the compile process
  G4Text dummytext = text;

  //----- 
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive( G4Text )\n";
#endif
  //----- Initialize DAWN IF NECESSARY
  FRBeginModeling();

  //----- send color
  /*
  const G4Color&	color     = GetTextColor (text) ;
  SendStrDouble3(	FR_COLOR_RGB    ,
			color.GetRed  (), 
			color.GetGreen(),
			color.GetBlue () );

  //----- Calc size 
  //MarkerSizeType size_type;
  G4double	fontsize  = GetMarkerDiameter( text , size_type );

  //----- Calc position
  const G4Point3D&	position       = text.GetPosition  () ;

  //----- offset
  G4double x_offset = text.GetXOffset();
  G4double y_offset = text.GetYOffset();

  //----- get string to be visualized and Calc its length
  const char* vis_text    = text.GetText();
  const int   STR_LENGTH  = strlen ( vis_text );

  //----- create buffer and copy the string there
  int   MAX_STR_LENGTH   =  COMMAND_BUF_SIZE - 100 ;
  if  ( MAX_STR_LENGTH <= 0 ) { 
    G4cerr << "ERROR (FukuiRenderer) : Not enough buffer size for data transferring." << G4endl;
    G4cerr << "                        G4Text Visualization is aborted" << G4endl;
    return ;
  }
  char*  buf  = new char [ (MAX_STR_LENGTH + 1) ] ; 
  if  ( MAX_STR_LENGTH >= STR_LENGTH  ) {
    strcpy  ( buf, vis_text ) ; 
  } else {
    strncpy ( buf, vis_text, MAX_STR_LENGTH ) ;
  }

  //----- select string command 
  char  text_command[32];
  switch (size_type) {
  case world:
    strcpy ( text_command, FR_MARK_TEXT_2D  );
    break;
  case screen:
  default:
    strcpy ( text_command, FR_MARK_TEXT_2DS );
    break;
  }

  //----- delete buffer
  delete [] buf ;
  */
} // G4GMocrenFileSceneHandler::AddPrimitive ( text )


//----- Add circle
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Circle& mark_circle )
{
  // to avoid a warning in the compile process
  G4Circle dummycircle = mark_circle;

  //----- 
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive( G4Circle )\n";
#endif
  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();


} // G4GMocrenFileSceneHandler::AddPrimitive ( mark_circle )


//----- Add square
void G4GMocrenFileSceneHandler::AddPrimitive (const G4Square& mark_square )
{
  // to avoid a warning in the compile process
  G4Square dummysquare = mark_square;

  //----- 
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive( G4Square )\n";
#endif
  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

} // G4GMocrenFileSceneHandler::AddPrimitive ( mark_square )


//----- Add polyhedron
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Polyhedron& polyhedron ) 
{
  //----- 
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddPrimitive( G4Polyhedron )\n";
#endif

  if (polyhedron.GetNoFacets() == 0) return;

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //---------- (3) Facet block
  for (int f = polyhedron.GetNoFacets(); f; f--){
    G4int notLastEdge;
    G4int index = -1; // initialization
    G4int edgeFlag = 1;
    G4int preedgeFlag = 1;
    G4int work[4], i = 0;
    do {
      preedgeFlag = edgeFlag;
      notLastEdge = polyhedron.GetNextVertexIndex(index, edgeFlag);
      work[i++] = index;
    }while (notLastEdge);
    switch (i){
    case 3:
      //SendStrInt3(FR_FACET, work[0], work[1], work[2] );
      break;
    case 4:
      //SendStrInt4(FR_FACET, work[0], work[1], work[2], work[3] );
      break;
    default:
      G4cerr <<
	"ERROR G4GMocrenFileSceneHandler::AddPrimitive(G4Polyhedron)\n";
      G4PhysicalVolumeModel* pPVModel =
        dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
      if (pPVModel) G4cerr <<
		      "Volume " << pPVModel->GetCurrentPV()->GetName() <<
		      ", Solid " << pPVModel->GetCurrentLV()->GetSolid()->GetName() <<
		      " (" << pPVModel->GetCurrentLV()->GetSolid()->GetEntityType();
      G4cerr <<
	"\nG4Polyhedron facet with " << i << " edges" << G4endl;	
    }
  }

} // G4GMocrenFileSceneHandler::AddPrimitive (polyhedron) 


//----- 
void G4GMocrenFileSceneHandler::FREndModeling ()
{
  G4VSceneHandler::EndModeling();

  //----- 		
#if defined DEBUG_FR_SCENE
  G4cerr << "***** FREndModeling (called)" << G4endl;
#endif
  if( FRIsInModeling() ) {

#if defined DEBUG_FR_SCENE
    G4cerr << "***** FREndModeling (started) " ; 
    G4cerr << "(/EndModeling, /DrawAll, /CloseDevice)" << G4endl;
#endif

    //----- End saving data to g4.gdd
    EndSavingGdd() ;

    //------ Reset flag 
    FRflag_in_modeling = false ;

  }

}


//----- 
void G4GMocrenFileSceneHandler::BeginPrimitives (const G4Transform3D& objectTransformation)
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** BeginPrimitives \n";
#endif

  FRBeginModeling();

  G4VSceneHandler::BeginPrimitives (objectTransformation);
  fpObjectTransformation = &objectTransformation;

}


//----- 
void G4GMocrenFileSceneHandler::EndPrimitives ()
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** EndPrimitives \n";
#endif
  G4VSceneHandler::EndPrimitives ();
}


//========== AddSolid() functions ==========//

//----- Add box
void G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( box )\n";
#endif

  if(GFDEBUG_DET > 0)
    G4cout << "G4GMocrenFileSceneHandler::AddSolid(const G4Box&)  : "
	   << box.GetName() << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();


  //--
  if(GFDEBUG_DET > 1) {
    G4cout << "-------" << G4endl;
    G4cout << "    " << box.GetName() << G4endl;
    G4Polyhedron * poly = box.CreatePolyhedron();
    poly->Transform(*fpObjectTransformation);
    //int nv = poly->GetNoVertices();
    G4Point3D v1, v2;
    G4int next;
    //while(1) { // next flag isn't functional.
    for(int i = 0; i < 12; i++) { // # of edges is 12.
      poly->GetNextEdge(v1, v2, next);
      if(next == 0) break;
      G4cout << "    (" << v1.x() << ", "
	     << v1.y() << ", "
	     << v1.z() << ") - ("
	     << v2.x() << ", "
	     << v2.y() << ", "
	     << v2.z() << ") [" << next << "]"
	     << G4endl;
    }
    delete poly;
  }


  const G4VModel* pv_model  = GetModel();
  if (!pv_model) { return ; } 
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) { return ; }


  //-- debug information
  if(GFDEBUG_DET > 0) {
    G4Material * mat = pPVModel->GetCurrentMaterial();
    G4String name = mat->GetName();
    G4double dens = mat->GetDensity()/(g/cm3);
    G4int copyNo = pPVModel->GetCurrentPV()->GetCopyNo();
    G4int depth = pPVModel->GetCurrentDepth();
    G4cout << "    copy no.: " << copyNo << G4endl;
    G4cout << "    depth   : " << depth << G4endl;
    G4cout << "    density : " << dens << " [g/cm3]" << G4endl;
    G4cout << "    location: " << pPVModel->GetCurrentPV()->GetObjectTranslation() << G4endl;
    G4cout << "    Multiplicity        : " << pPVModel->GetCurrentPV()->GetMultiplicity() << G4endl;
    G4cout << "    Is replicated?      : " << pPVModel->GetCurrentPV()->IsReplicated() << G4endl;
    G4cout << "    Is parameterised?   : " << pPVModel->GetCurrentPV()->IsParameterised() << G4endl;
    G4cout << "    top phys. vol. name : " << pPVModel->GetTopPhysicalVolume()->GetName() << G4endl;
  }

  //-- check parameterised volume
  G4String volName = fMessenger.getVolumeName();
  if(box.GetName() == volName) {
  //G4String paraName = "DICOMNest";
  //if(box.GetName() == paraName) {
    
    fVolumeTrans3D = *fpObjectTransformation;
    // coordination system correction for gMocren
    G4ThreeVector raxis(1., 0., 0.), dummy(0.,0.,0.); 
    G4RotationMatrix rot(raxis, M_PI*rad);
    G4Transform3D trot(rot, dummy);
    if(GFDEBUG_DET) {
      G4ThreeVector trans = fVolumeTrans3D.getTranslation();
      G4RotationMatrix rot = fVolumeTrans3D.getRotation().inverse();
      G4cout << "fVolumeTrans3D: " << trans << G4endl << rot << G4endl;
    }
    fVolumeTrans3D = fVolumeTrans3D*trot;
    G4cout << " Parameterised volume : " << box.GetName() << G4endl;



    //
    G4VPhysicalVolume * pv[3] = {0,0,0};
    pv[0] = pPVModel->GetCurrentPV()->GetLogicalVolume()->GetDaughter(0);
    if(!pv[0]) {
      G4Exception("Error[gMocrenFileSceneHandler]: Unexpected volume.");
    }
    G4int dirAxis[3] = {-1,-1,-1};
    G4int nDaughters[3] = {0,0,0};

    EAxis axis; G4int nReplicas; G4double width; G4double offset; G4bool consuming;  
    pv[0]->GetReplicationData(axis, nReplicas, width, offset, consuming);
    nDaughters[0] = nReplicas;
    switch(axis) {
    case kXAxis: dirAxis[0] = 0; break;
    case kYAxis: dirAxis[0] = 1; break;
    case kZAxis: dirAxis[0] = 2; break;
    default: G4Exception("Error.");
    }
    fNestedVolumeNames.push_back(pv[0]->GetName());
    if(GFDEBUG_DET) 
      G4cout << "        daughter name :  " << pv[0]->GetName()
	     << "   # : " << nDaughters[0] << G4endl;

    //
    pv[1] = pv[0]->GetLogicalVolume()->GetDaughter(0);
    if(pv[1]) {
      pv[1]->GetReplicationData(axis, nReplicas, width, offset, consuming);
      nDaughters[1] = nReplicas;
      switch(axis) {
      case kXAxis: dirAxis[1] = 0; break;
      case kYAxis: dirAxis[1] = 1; break;
      case kZAxis: dirAxis[1] = 2; break;
      default: G4Exception("Error.");
      }
      fNestedVolumeNames.push_back(pv[1]->GetName());
      if(GFDEBUG_DET) 
	G4cout << "        sub-daughter name :  " << pv[1]->GetName()
	      << "   # : " << nDaughters[1]<< G4endl;

      //
      pv[2] = pv[1]->GetLogicalVolume()->GetDaughter(0);
      if(pv[2]) {
	nDaughters[2] = pv[2]->GetMultiplicity();
	fNestedVolumeNames.push_back(pv[2]->GetName());
	if(GFDEBUG_DET) 
	  G4cout << "        sub-sub-daughter name :  " << pv[2]->GetName()
		 << "   # : " << nDaughters[2] << G4endl;

	if(nDaughters[2] > 1) {
	  G4VNestedParameterisation * nestPara
	    = dynamic_cast<G4VNestedParameterisation*>(pv[2]->GetParameterisation());
	  if(!nestPara) G4Exception("Error[gMocrenFileSceneHandler]: None nested parameterisation");
	  nestPara->ComputeTransformation(0, pv[2]);
	  G4ThreeVector trans0 = pv[2]->GetObjectTranslation();
	  nestPara->ComputeTransformation(1, pv[2]);
	  G4ThreeVector trans1 = pv[2]->GetObjectTranslation();
	  G4ThreeVector diff(trans0 - trans1);
	  if(GFDEBUG_DET) 
	    G4cout << trans0 << " - " << trans1 << " - " << diff << G4endl;

	  if(diff.x() != 0.) dirAxis[2] = 0;
	  else if(diff.y() != 0.) dirAxis[2] = 1;
	  else if(diff.z() != 0.) dirAxis[2] = 2;
	  else G4Exception("Error[gMocrenFileSceneHandler]: Unexpected nested parameterisation");
	}
      }
    }

    // get densities
    G4VNestedParameterisation * nestPara
      = dynamic_cast<G4VNestedParameterisation*>(pv[2]->GetParameterisation());
    if(nestPara) {
      G4double prexyz[3] = {0.,0.,0.}, xyz[3] = {0.,0.,0.};
      for(int n0 = 0; n0 < nDaughters[0]; n0++) {
	for(int n1 = 0; n1 < nDaughters[1]; n1++) {
	  for(int n2 = 0; n2 < nDaughters[2]; n2++) {
		  
	    G4GMocrenTouchable * touch = new G4GMocrenTouchable(n1, n0);
	    if(GFDEBUG_DET) 
	      G4cout << "   retrieve volume : copy # : " << n0
		     << ", " << n1 << ", " << n2 << G4endl;
	    G4Material * mat = nestPara->ComputeMaterial(pv[2], n2, touch);
	    delete touch;
	    G4double dens = mat->GetDensity()/(g/cm3);
 
	    if(GFDEBUG_DET) 
	      G4cout << "           density :" << dens << " [g/cm3]" << G4endl;
	    G4Box tbox(box);
	    nestPara->ComputeDimensions(tbox, n2, pv[2]);
	    xyz[0] = tbox.GetXHalfLength()/mm;
	    xyz[1] = tbox.GetYHalfLength()/mm;
	    xyz[2] = tbox.GetZHalfLength()/mm;
	    if(n0 != 0 || n1 != 0 || n2 != 0) {
	      for(int i = 0; i < 3; i++) {
		if(xyz[i] != prexyz[i]) G4Exception("Error[gMocrenFileSceneHandler]: Unsupported parameterisation.");
	      }
	    }
	    if(GFDEBUG_DET) 
	      G4cout << "              size : " << tbox.GetXHalfLength()/mm << " x "
		     << tbox.GetYHalfLength()/mm << " x "
		     << tbox.GetZHalfLength()/mm << " [mm3]" << G4endl;

	    G4int idx[3];
	    idx[dirAxis[0]] = n0;
	    idx[dirAxis[1]] = n1;
	    idx[dirAxis[2]] = n2;
	    Index3D i3d(idx[0],idx[1],idx[2]);
	    fNestedModality[i3d] = dens;
	    if(GFDEBUG_DET) 
	      G4cout << " index: " << idx[0] << ", " << idx[1] << ", " << idx[2]
		     << "  density: " << dens << G4endl;

	    for(int i = 0; i < 3; i++) prexyz[i] = xyz[i];
	  }
        }
      }  

      fVolumeSize.set(box.GetXHalfLength()*2/mm,
		      box.GetYHalfLength()*2/mm,
		      box.GetZHalfLength()*2/mm);
      // mesh size
      if(!fbSetModalityVoxelSize) {
	float spacing[3] = {2*xyz[0], 2*xyz[1], 2*xyz[2]};
	fgMocrenIO->setVoxelSpacing(spacing);
	fVoxelDimension.set(spacing[0], spacing[1], spacing[2]);
	fbSetModalityVoxelSize = true;
      }
    } else {
      if(GFDEBUG_DET) 
        G4cout << pv[2]->GetName() << G4endl;
      G4Exception("Error[gMocrenFileSceneHandler]: none nested parameterization");
    }
  }


  //-- add detectors
  G4bool bAddDet = true;
  if(!fMessenger.getDrawVolumeGrid()) {

    if(volName == box.GetName()) {
      bAddDet = false;
    }

    std::vector<G4String>::iterator itr = fNestedVolumeNames.begin();
    for(; itr != fNestedVolumeNames.end(); itr++) {
      if(*itr == box.GetName())  {
	bAddDet = false;
	break;
      }
    }

  }
  if(bAddDet) AddDetector(box);
    


  //-- debug
  if(GFDEBUG_DET > 1) {
    if(pPVModel->GetCurrentPV()->IsParameterised()) {
      G4VPVParameterisation * para = pPVModel->GetCurrentPV()->GetParameterisation();
      G4cout << " Is nested parameterisation? : " << para->IsNested() << G4endl;


      G4int npvp = pPVModel->GetDrawnPVPath().size();
      G4cout << "     physical volume node id : "
	     << "size: " << npvp << ", PV name: ";
      for(int i = 0; i < npvp; i++) {
	G4cout << pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->GetName()
	       << " [param:"
	       << pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->IsParameterised()
	       << ",rep:"
	       << pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->IsReplicated();
	if(pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->GetParameterisation()) {
	  G4cout << ",nest:"
		 << pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->GetParameterisation()->IsNested();
	}
	G4cout << ",copyno:"
	       << pPVModel->GetDrawnPVPath()[i].GetPhysicalVolume()->GetCopyNo();
	G4cout << "] - ";
      }
      G4cout << G4endl;


      EAxis axis; G4int nReplicas; G4double width; G4double offset; G4bool consuming;  
      pPVModel->GetCurrentPV()->GetReplicationData(axis, nReplicas, width, offset, consuming);
      G4cout << "     # replicas : " << nReplicas << G4endl;
      G4double pareDims[3];
      G4Box * pbox = dynamic_cast<G4Box *>(pPVModel->GetDrawnPVPath()[npvp-2].GetPhysicalVolume()->GetLogicalVolume()->GetSolid());
      if(pbox) {
	pareDims[0] = 2.*pbox->GetXHalfLength()*mm;
	pareDims[1] = 2.*pbox->GetYHalfLength()*mm;
	pareDims[2] = 2.*pbox->GetZHalfLength()*mm;
	G4cout << "     mother size ["
	       << pPVModel->GetDrawnPVPath()[npvp-2].GetPhysicalVolume()->GetName()
	       << "] : "
	       << pareDims[0] << " x "
	       << pareDims[1] << " x "
	       << pareDims[2] << " [mm3]"
	       << G4endl;
      }
      G4double paraDims[3];
      G4Box * boxP = dynamic_cast<G4Box *>(pPVModel->GetDrawnPVPath()[npvp-1].GetPhysicalVolume()->GetLogicalVolume()->GetSolid());
      if(boxP) {
	paraDims[0] = 2.*boxP->GetXHalfLength()*mm;
	paraDims[1] = 2.*boxP->GetYHalfLength()*mm;
	paraDims[2] = 2.*boxP->GetZHalfLength()*mm;
	G4cout << "     parameterised volume? ["
	       << pPVModel->GetDrawnPVPath()[npvp-1].GetPhysicalVolume()->GetName()
	       << "] : "
	       << paraDims[0] << " x "
	       << paraDims[1] << " x "
	       << paraDims[2] << " [mm3]  : "
	       << G4int(pareDims[0]/paraDims[0]) << " x "
	       << G4int(pareDims[1]/paraDims[1]) << " x "
	       << G4int(pareDims[2]/paraDims[2]) << G4endl;
      } else {
	G4cout << pPVModel->GetDrawnPVPath()[npvp-2].GetPhysicalVolume()->GetName()
	       << " isn't a G4Box." << G4endl;
      }
    }
  }




} // void G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )


//----- Add tubes
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Tubs& tubes )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( tubes )\n";
#endif
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(tubes);


  //AK test
  if(GFDEBUG_DET > 0) {
    G4cout << "-------" << G4endl;
    G4cout << "    " << tubes.GetName() << G4endl;
    G4Polyhedron * poly = tubes.CreatePolyhedron();
    int nv = poly->GetNoVertices();
    for(int i = 0; i < nv; i++) {
      G4cout << "    (" << poly->GetVertex(i).x() << ", "
	     << poly->GetVertex(i).y() << ", "
	     << poly->GetVertex(i).z() << ")" << G4endl;
    }
    delete poly;
  }

  const G4VModel* pv_model  = GetModel();
  if (!pv_model) { return ; } 
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) { return ; }
  G4Material * mat = pPVModel->GetCurrentMaterial();
  G4String name = mat->GetName();
  //G4double dens = mat->GetDensity() /(g/cm3);
  //SendStrDouble("/DENSITY", dens);
} // void G4GMocrenFileSceneHandler::AddSolid( const G4Tubs& )



//----- Add cons
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Cons& cons )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( cons )\n";
#endif
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(cons);

}// G4GMocrenFileSceneHandler::AddSolid( cons )


//----- Add trd
void G4GMocrenFileSceneHandler::AddSolid ( const G4Trd& trd )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( trd )\n";
#endif

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(trd);

} // G4GMocrenFileSceneHandler::AddSolid ( trd )


//----- Add sphere
void G4GMocrenFileSceneHandler::AddSolid ( const G4Sphere& sphere )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( sphere )\n";
#endif
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(sphere);

} // G4GMocrenFileSceneHandler::AddSolid ( sphere )


//----- Add para
void G4GMocrenFileSceneHandler::AddSolid (const G4Para& para)
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( para )\n";
#endif

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(para);

} // G4GMocrenFileSceneHandler::AddSolid ( para )


//----- Add trap
void G4GMocrenFileSceneHandler::AddSolid (const G4Trap& trap)
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( trap )\n";
#endif

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(trap);

} // G4GMocrenFileSceneHandler::AddSolid (const G4Trap& trap)


//----- Add torus
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Torus& torus )
{
#if defined DEBUG_FR_SCENE
  G4cerr << "***** AddSolid ( torus )\n";
#endif
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(torus);

} // void G4GMocrenFileSceneHandler::AddSolid( const G4Torus& )



//----- Add a shape which is not treated above
void G4GMocrenFileSceneHandler::AddSolid ( const G4VSolid& solid  )
{
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize Fukui Renderer IF NECESSARY
  FRBeginModeling();

  //
  AddDetector(solid);

  //----- Send a primitive
  G4VSceneHandler::AddSolid( solid ) ; 

} //G4GMocrenFileSceneHandler::AddSolid ( const G4VSolid& ) 

#include "G4TrajectoriesModel.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"

//----- Add a trajectory
void G4GMocrenFileSceneHandler::AddCompound(const G4VTrajectory & traj) {

  fbModelingTrajectory = true;

  G4VSceneHandler::AddCompound(traj);

  if(GFDEBUG_TRK) {
    std::cout << " ::AddCompound(const G4VTrajectory&) >>>>>>>>> " << std::endl;

    G4TrajectoriesModel * pTrModel = dynamic_cast<G4TrajectoriesModel*>(fpModel);
    if (!pTrModel) { 
      G4Exception 
	("G4VSceneHandler::AddCompound(const G4VTrajectory&): Not a G4TrajectoriesModel.");
    } else {
      traj.DrawTrajectory(pTrModel->GetDrawingMode());

      const G4VTrajectory * trj = pTrModel->GetCurrentTrajectory();
      G4cout << "------ track" << G4endl;
      G4cout << "    name:     " << trj->GetParticleName() << G4endl;
      G4cout << "    id:       " << trj->GetTrackID() << G4endl;
      G4cout << "    charge:   " << trj->GetCharge() << G4endl;
      G4cout << "    momentum: " << trj->GetInitialMomentum() << G4endl;
      
      int nPnt = trj->GetPointEntries();
      G4cout << "    point:    ";
      for(int i = 0; i < nPnt; i++) {
	G4cout << trj->GetPoint(i)->GetPosition() << ", ";
      }
      G4cout << G4endl;


    }
  }

  fbModelingTrajectory = false;
}

#include <vector>
#include "G4VHit.hh"
#include "G4AttValue.hh"
//----- Add a hit 
void G4GMocrenFileSceneHandler::AddCompound( const G4VHit & hit) {
  if(GFDEBUG_HIT) G4cout << " ::AddCompound(const G4VHit&) >>>>>>>>> " << G4endl;

  G4VSceneHandler::AddCompound(hit);

  /*
    const std::map<G4String, G4AttDef> * map = hit.GetAttDefs();
    if(!map) return;
    std::map<G4String, G4AttDef>::const_iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
    G4cout << itr->first << " : " << itr->second.GetName()
    << " , " << itr->second.GetDesc() << G4endl;
    }
  */

  std::vector<G4String> hitNames = fMessenger.getHitNames();
  if(GFDEBUG_HIT) {
    std::vector<G4String>::iterator itr = hitNames.begin();
    for(; itr != hitNames.end(); itr++) 
      G4cout << "  hit name : " << *itr << G4endl;
  }
  
  std::vector<G4AttValue> * attval = hit.CreateAttValues();
  if(!attval) {G4cout << "0 empty " << (unsigned long)attval << G4endl;}
  else {

    G4bool bid[3] = {false, false, false};
    Index3D id;

    G4int nhitname = (G4int)hitNames.size();

    std::vector<G4AttValue>::iterator itr;
    // First, get IDs
    for(itr = attval->begin(); itr != attval->end(); itr++) {
      std::string stmp = itr->GetValue();
      std::istringstream sval(stmp.c_str());

      if(itr->GetName() == G4String("XID")) {
	sval >>	id.x;
	bid[0] = true;
	continue;
      }
      if(itr->GetName() == G4String("YID")) {
	sval >>	id.y;
	bid[1] = true;
	continue;
      }
      if(itr->GetName() == G4String("ZID")) {
	sval >>	id.z;
	bid[2] = true;
	continue;
      }
    }

    if(bid[0] && bid[1] && bid[2]) {

      if(GFDEBUG_HIT)
	G4cout << " Hit : index(" << id.x << ", " << id.y << ", "
	       << id.z << ")" << G4endl;

      // Get attributes
      for(itr = attval->begin(); itr != attval->end(); itr++) {
	for(int i = 0; i < nhitname; i++) {
	  if(itr->GetName() == hitNames[i]) {

	    std::string stmp = itr->GetValue();
	    std::istringstream sval(stmp.c_str());
	    G4double value;
	    G4String unit;
	    sval >> value >> unit;

	    std::map<G4String, std::map<Index3D, G4double> >::iterator fNestedHitsListItr;
	    fNestedHitsListItr = fNestedHitsList.find(hitNames[i]);
	    if(fNestedHitsListItr != fNestedHitsList.end()) {
	      //fTempNestedHits = &fNestedHitsListItr->second;
	      //(*fTempNestedHits)[id] = value;
	      fNestedHitsListItr->second[id] = value;
	    } else {
	      std::map<Index3D, G4double> hits;
	      hits[id] = value;
	      fNestedHitsList[hitNames[i]] = hits;
	    }

	    
	    if(GFDEBUG_HIT)
	      G4cout << "     : " << hitNames[i] << " -> " << value
		     << " [" << unit << "]" << G4endl;
	  }
	}
      }
    } else {
      G4Exception("Error in G4GMocrenFileSceneHandler::AddCompound(const G4VHit &)");
    }

    delete attval;
  }

}

void G4GMocrenFileSceneHandler::AddCompound(const G4THitsMap<G4double> & hits) {
  if(GFDEBUG_HIT)
    G4cout << " ::AddCompound(const std::map<G4int, G4double*> &) >>>>>>>>> " << G4endl;

  if(GFDEBUG_HIT) {
    G4String meshname = ((G4VHitsCollection)hits).GetSDname();
    G4String scorername = ((G4VHitsCollection)hits).GetName();
    G4cout << "       >>>>> " << meshname << " : " << scorername  << G4endl;

    std::vector<G4String> hitScorerNames = fMessenger.getHitScorerNames();
    G4int nhitname = (G4int)hitScorerNames.size();
    for(int i = 0; i < nhitname; i++)
      if(scorername == hitScorerNames[i]) 
	G4cout << "       !!!! Hit scorer !!!!" << G4endl;

    std::map<G4int, G4double*> * map = hits.GetMap();
    std::map<G4int, G4double*>::const_iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      G4cout << "[" << itr->first << "] " << *itr->second << ", ";
    }
    G4cout << G4endl;
  }

}

//----- 
G4bool G4GMocrenFileSceneHandler::IsVisible()
{
  //----- 
  G4bool  visibility  = true ;

  //----- 
  const G4VisAttributes* pVisAttribs =
    fpViewer->GetApplicableVisAttributes( fpVisAttribs );

  //----- 
  if( ( getenv( FR_ENV_CULL_INVISIBLE_OBJECTS ) != NULL      )   && \
      ( strcmp( getenv( FR_ENV_CULL_INVISIBLE_OBJECTS ),"0"  ) ) && \
      ( pVisAttribs )                                             ) 
    {
      visibility = pVisAttribs->IsVisible();
    } 

  //----- 
  return visibility ;

} // G4GMocrenFileSceneHandler::IsVisible()

//----- 
void G4GMocrenFileSceneHandler::ClearTransientStore() 
{
  G4VSceneHandler::ClearTransientStore ();
  // This is typically called after an update and before drawing hits
  // of the next event.  To simulate the clearing of "transients"
  // (hits, etc.) the detector is redrawn...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

//----- 
void G4GMocrenFileSceneHandler::AddDetector(const G4VSolid & solid) {

  Detector detector;

  // detector name
  detector.name = solid.GetName();
  if(GFDEBUG_DET > 1)
    G4cout << "0 Detector name : " << detector.name << G4endl;
 
  const G4VModel* pv_model  = GetModel();
  if (!pv_model) { return ; } 
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) { return ; }

  // edge points of the detector
  std::vector<float *> dedges;
  G4Polyhedron * poly = solid.CreatePolyhedron();
  detector.polyhedron = poly;
  detector.transform3D = *fpObjectTransformation;

  // retrieve color
  unsigned char uccolor[3] = {30, 30, 30};
  if(pPVModel->GetCurrentLV()->GetVisAttributes()) {
    G4Color color = pPVModel->GetCurrentLV()->GetVisAttributes()->GetColor();
    uccolor[0] = (unsigned char)(color.GetRed()*255);
    uccolor[1] = (unsigned char)(color.GetGreen()*255);
    uccolor[2] = (unsigned char)(color.GetBlue()*255);
    //if(uccolor[0] < 2 && uccolor[1] < 2 && uccolor[2] < 2)
    //uccolor[0] = uccolor[1] = uccolor[2] = 30; // dark grey
  }
  for(int i = 0; i < 3; i++) detector.color[i] = uccolor[i];
  //
  fDetectors.push_back(detector);

  if(GFDEBUG_DET > 1) {
    G4cout << "0     color:   (" << (int)uccolor[0] << ", "
	   << (int)uccolor[1] << ", " << (int)uccolor[2] << ")"
	   << G4endl;
  }

}

//----- 
void G4GMocrenFileSceneHandler::ExtractDetector() {

  std::vector<Detector>::iterator itr = fDetectors.begin();

  for(; itr != fDetectors.end(); itr++) {

    // detector name
    G4String detname = itr->name;
    if(GFDEBUG_DET > 0)
      G4cout << "Detector name : " << detname << G4endl;

    // edge points of the detector
    std::vector<float *> dedges;
    G4Polyhedron * poly = itr->polyhedron;
    poly->Transform(itr->transform3D);
    G4Transform3D invVolTrans = fVolumeTrans3D.inverse();
    poly->Transform(invVolTrans);

    G4Point3D v1, v2;
    G4bool bnext = true;
    G4int next;
    G4int nedges = 0;
    //
    while(bnext) {
      if(!(poly->GetNextEdge(v1, v2, next))) bnext =false;
      float * edge = new float[6];
      edge[0] = v1.x()/mm;
      edge[1] = v1.y()/mm;
      edge[2] = v1.z()/mm;
      edge[3] = v2.x()/mm;
      edge[4] = v2.y()/mm;
      edge[5] = v2.z()/mm;
      dedges.push_back(edge);
      nedges++;
    }
    //delete poly;
    // detector color
    unsigned char uccolor[3] = {itr->color[0],
				itr->color[1],
				itr->color[2]};
    //
    fgMocrenIO->addDetector(detname, dedges, uccolor);
    for(int i = 0; i < nedges; i++) { // # of edges is 12.
      delete [] dedges[i];
    }
    dedges.clear(); 

    if(GFDEBUG_DET > 0) {
      G4cout << "    color:   (" << (int)uccolor[0] << ", "
	     << (int)uccolor[1] << ", " << (int)uccolor[2] << ")"
	     << G4endl;
    }
  }
}

G4GMocrenFileSceneHandler::Detector::Detector()
  : polyhedron(0) {
  color[0] = color[1] = color[2] = 255;
}
G4GMocrenFileSceneHandler::Detector::~Detector() {
  if(!polyhedron) delete polyhedron;
}
void G4GMocrenFileSceneHandler::Detector::clear() {
  name.clear();
  if(!polyhedron) delete polyhedron;
  color[0] = color[1] = color[2] = 255;
  transform3D = G4Transform3D::Identity;
}

G4GMocrenFileSceneHandler::Index3D::Index3D(G4int _x, G4int _y, G4int _z) 
  : x(_x), y(_y), z(_z) {
  ;
}

G4GMocrenFileSceneHandler::Index3D::Index3D()
  : x(0), y(0), z(0) {
}
G4bool G4GMocrenFileSceneHandler::Index3D::operator < (const Index3D & _right) const {
  if(z < _right.z) {
     return true;
  } else if(z == _right.z) {
    if(y < _right.y) return true;
    else if(y == _right.y) 
      if(x < _right.x) return true;
  } 
  return false;
}
G4bool G4GMocrenFileSceneHandler::Index3D::operator == (const Index3D & _right) const {
  if(z == _right.z && y == _right.y && x == _right.x) return true;
  return false;
}



//////////////////////
// static variables //
//////////////////////

//----- static variables
G4int G4GMocrenFileSceneHandler::fSceneIdCount = 0; 

