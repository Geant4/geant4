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
// $Id: G4GMocrenFileSceneHandler.cc 104291 2017-05-23 13:25:51Z gcosmo $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//           Sep. 22, 2009  Akinori Kimura : modify and fix code to support 
//                                          PrimitiveScorers and clean up
//
// GMocrenFile scene handler


//----- header files
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <sstream>
#include <iomanip>

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisManager.hh"

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
#include "G4TrajectoriesModel.hh"
#include "G4VTrajectoryModel.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "G4HitsModel.hh"
#include "G4GMocrenMessenger.hh"
//#include "G4PSHitsModel.hh"
#include "G4GMocrenIO.hh"
#include "G4VNestedParameterisation.hh"
#include "G4GMocrenTouchable.hh"
#include "G4GMocrenFileCTtoDensityMap.hh"
#include "G4PhantomParameterisation.hh"
#include "G4PhysicalVolumeSearchScene.hh"

#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"

//----- constants
const char  GDD_FILE_HEADER      [] = "g4_";
const char  DEFAULT_GDD_FILE_NAME[] = "g4_00.gdd";

const G4int FR_MAX_FILE_NUM = 100 ;
const G4int MAX_NUM_TRAJECTORIES = 100000;

//-- for a debugging
const G4bool GFDEBUG = false;
const G4bool GFDEBUG_TRK = false;//true;
const G4bool GFDEBUG_HIT = false;//true;
const G4bool GFDEBUG_DIGI = false;//true;
const G4int GFDEBUG_DET = 0; // 0: false 

//////////////////////
// static variables //
//////////////////////

//----- static variables
G4int G4GMocrenFileSceneHandler::kSceneIdCount = 0; 

///////////////////////////
// Driver-dependent part //
///////////////////////////


//----- G4GMocrenFileSceneHandler, constructor
G4GMocrenFileSceneHandler::G4GMocrenFileSceneHandler(G4GMocrenFile& system,
						     G4GMocrenMessenger & messenger,
						     const G4String& name)
  : G4VSceneHandler(system, kSceneIdCount++, name),
    kSystem(system),
    kMessenger(messenger),
    kgMocrenIO(new G4GMocrenIO()),
    kbSetModalityVoxelSize(false),
    kbModelingTrajectory(false),
//    kGddDest(0),
    kFlagInModeling(false),
    kFlagSaving_g4_gdd(false),
    kFlagParameterization(0),
    kFlagProcessedInteractiveScorer(false) {

  // g4.gdd filename and its directory
  if(std::getenv("G4GMocrenFile_DEST_DIR") == NULL) {
    kGddDestDir[0] = '\0';
    //std::strcpy(kGddDestDir , "");                    // output dir
    //std::strcpy(kGddFileName, DEFAULT_GDD_FILE_NAME); // filename
    std::strncpy(kGddFileName, DEFAULT_GDD_FILE_NAME,
		 std::strlen(DEFAULT_GDD_FILE_NAME)+1);   // filename
  } else {
    const char * env = std::getenv("G4GMocrenFile_DEST_DIR");
    int len = std::strlen(env);
    if(len > 256) {
      G4Exception("G4GMocrenFileSceneHandler::G4GMocrenFileSceneHandler(*)",
                  "gMocren1000", FatalException,
                  "Invalid length of string set in G4GMocrenFile_DEST_DIR");
    }
    std::strncpy(kGddDestDir, env, len+1);  // output dir
    std::strncpy(kGddFileName, DEFAULT_GDD_FILE_NAME,
                 std::strlen(DEFAULT_GDD_FILE_NAME)+1);  // filename
  }
		
  // maximum number of g4.gdd files in the dest directory
  kMaxFileNum = FR_MAX_FILE_NUM ; // initialization
  if ( std::getenv( "G4GMocrenFile_MAX_FILE_NUM" ) != NULL ) {	
    char * pcFileNum = getenv("G4GMocrenFile_MAX_FILE_NUM");
    char c10FileNum[10];
    std::strncpy(c10FileNum, pcFileNum, 10);
    kMaxFileNum = std::atoi(c10FileNum);

  } else {
    kMaxFileNum = FR_MAX_FILE_NUM ;
  }
  if( kMaxFileNum < 1 ) { kMaxFileNum = 1 ; }

  InitializeParameters();

} 


//----- G4GMocrenFileSceneHandler, destructor
G4GMocrenFileSceneHandler::~G4GMocrenFileSceneHandler () 
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "***** ~G4GMocrenFileSceneHandler" << G4endl;

  if(kGddDest) {
    //----- End of modeling
    // close g4.gdd
    GFEndModeling();
  }
  if(kgMocrenIO != NULL) delete kgMocrenIO;

}

//----- initialize all parameters
void G4GMocrenFileSceneHandler::InitializeParameters() {

  kbSetModalityVoxelSize = false;

  for(G4int i = 0; i < 3; i++) {
    kModalitySize[i] = 0;
    kNestedVolumeDimension[i] = 0;
    kNestedVolumeDirAxis[i] = -1;
  }

  // delete kgMocrenIO;

}

//-----
void G4GMocrenFileSceneHandler::SetGddFileName() 
{
  // g4_00.gdd, g4_01.gdd, ..., g4_MAX_FILE_INDEX.gdd
  const G4int MAX_FILE_INDEX = kMaxFileNum - 1 ;

  // dest directory (null if no environmental variables is set)
  std::strncpy(kGddFileName, kGddDestDir, std::strlen(kGddDestDir)+1);

  // create full path name (default)
  std::strncat ( kGddFileName, DEFAULT_GDD_FILE_NAME, std::strlen(DEFAULT_GDD_FILE_NAME));

  // Automatic updation of file names
  static G4int currentNumber = 0;
  for( G4int i = currentNumber ; i < kMaxFileNum ; i++) { 

    // Message in the final execution
    if( i == MAX_FILE_INDEX ) 
      {
	if (G4VisManager::GetVerbosity() >= G4VisManager::warnings) {
	  G4cout << "==========================================="   << G4endl; 
	  G4cout << "WARNING MESSAGE from GMocrenFile driver:   "   << G4endl;
	  G4cout << "  This file name is the final one in the   "   << G4endl;
	  G4cout << "  automatic updation of the output file name." << G4endl; 
	  G4cout << "  You may overwrite existing files, i.e.   "   << G4endl; 
	  G4cout << "  g4_XX.gdd."   << G4endl;
	  G4cout << "==========================================="   << G4endl; 
	}
      }

    // re-determine file name as G4GMocrenFile_DEST_DIR/g4_XX.gdd 
    std::ostringstream filename;
    filename
    << kGddDestDir << GDD_FILE_HEADER
    << std::setw(2) << std::setfill('0') << i << ".wrl";
    strncpy(kGddFileName,filename.str().c_str(),sizeof(kGddFileName));

    // check validity of the file name
    std::ifstream fin(kGddFileName); 
    if(GFDEBUG)
      G4cout << "FILEOPEN: " << i << " : " << kGddFileName << fin.fail()
	     << G4endl;
    if(!fin) { 
      // new file	
      fin.close();
      currentNumber = i+1;
      break; 
    } else { 
      // already exists (try next) 
      fin.close(); 
    } 

  } // for 

  G4cout << "======================================================================" << G4endl; 
  G4cout << "Output file: " << kGddFileName                          << G4endl; 
  G4cout << "Destination directory (current dir if NULL): " << kGddDestDir << G4endl; 
  G4cout << "Maximum number of files in the destination directory: " << kMaxFileNum << G4endl; 
  G4cout << "Note:" << G4endl; 
  G4cout << "  * The maximum number is customizable as:           " << G4endl;
  G4cout << "      % setenv  G4GMocrenFile_MAX_FILE_NUM  number " << G4endl;        
  G4cout << "  * The destination directory is customizable as:" << G4endl;
  G4cout << "      % setenv  G4GMocrenFile_DEST_DIR  dir_name/  " << G4endl;        
  G4cout << "     ** Do not forget \"/\" at the end of the dir_name, e.g. \"./tmp/\"." << G4endl;              
  //G4cout << "        dir_name, e.g. \"./tmp/\"."                 << G4endl;              
  G4cout << G4endl;
  G4cout << "Maximum number of trajectories is set to " << MAX_NUM_TRAJECTORIES << "."<< G4endl;
  G4cout << "======================================================================" << G4endl; 

} // G4GMocrenFileSceneHandler::SetGddFileName()


//-----
void	G4GMocrenFileSceneHandler::BeginSavingGdd( void )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "***** BeginSavingGdd (called)" << G4endl;

  if( !IsSavingGdd() ) {

    if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations) {
      G4cout << "*****                   (started) " ;
      G4cout << "(open g4.gdd, ##)"  << G4endl;
    }

    SetGddFileName() ; // result set to kGddFileName
    kFlagSaving_g4_gdd = true; 


    G4GMocrenFileCTtoDensityMap ctdens;
    short minmax[2];
    minmax[0] = ctdens.GetMinCT();
    minmax[1] = ctdens.GetMaxCT();
    kgMocrenIO->setModalityImageMinMax(minmax);
    std::vector<G4float> map;
    G4float dens;
    for(G4int i = minmax[0]; i <= minmax[1]; i++) {
      dens = ctdens.GetDensity(i);
      map.push_back(dens);
    }
    kgMocrenIO->setModalityImageDensityMap(map);

    /*
    G4String fname = "modality-map.dat";
    std::ifstream ifile(fname);
    if(ifile) {
      short minmax[2];
      ifile >> minmax[0] >> minmax[1];
      kgMocrenIO->setModalityImageMinMax(minmax);
      std::vector<G4float> map;
      G4float dens;
      for(G4int i = minmax[0]; i <= minmax[1]; i++) {
	ifile >> dens;
	map.push_back(dens);
      }
      kgMocrenIO->setModalityImageDensityMap(map);
      
    } else {
      G4cout << "cann't open the file : " << fname << G4endl;
    }
    */

    // mesh size
    //kMessenger.getNoVoxels(kModalitySize[0], kModalitySize[1], kModalitySize[2]);
    //kgMocrenIO->setModalityImageSize(kModalitySize);
    
    // initializations
    //kgMocrenIO->clearModalityImage();
    kgMocrenIO->clearDoseDistAll();
    kgMocrenIO->clearROIAll();
    kgMocrenIO->clearTracks();
    kgMocrenIO->clearDetector();
    std::vector<Detector>::iterator itr = kDetectors.begin();
    for(; itr != kDetectors.end(); itr++) {
      itr->clear();
    }
    kDetectors.clear();
    
    kNestedHitsList.clear();
    kNestedVolumeNames.clear();
      
  }
}

void	G4GMocrenFileSceneHandler::EndSavingGdd  ( void ) 
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** EndSavingGdd (called)" << G4endl;

  if(IsSavingGdd()) {
    if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "*****                 (started) (close "
	     << kGddFileName << ")" << G4endl;

    if(kGddDest) kGddDest.close();
    kFlagSaving_g4_gdd = false; 

    std::map<Index3D, G4float>::iterator itr = kNestedModality.begin();
    G4int xmax=0, ymax=0, zmax=0;
    for(; itr != kNestedModality.end(); itr++) {
      if(itr->first.x > xmax) xmax = itr->first.x;
      if(itr->first.y > ymax) ymax = itr->first.y;
      if(itr->first.z > zmax) zmax = itr->first.z;
    }
    // mesh size
    kModalitySize[0] = xmax+1;
    kModalitySize[1] = ymax+1;
    kModalitySize[2] = zmax+1;
    kgMocrenIO->setModalityImageSize(kModalitySize);
    if(GFDEBUG) G4cout << "gMocren-file driver : modality size : "
		       << kModalitySize[0] << " x "
		       << kModalitySize[1] << " x "
		       << kModalitySize[2] << G4endl;

    G4int nxy = kModalitySize[0]*kModalitySize[1];
    //std::map<G4int, G4float>::iterator itr;
    for(G4int z = 0; z < kModalitySize[2]; z++) {
      short * modality = new short[nxy];
      for(G4int y = 0; y < kModalitySize[1]; y++) {
	for(G4int x = 0; x < kModalitySize[0]; x++) {
	  //for(G4int x = kModalitySize[0]-1; x >= 0 ; x--) {
	  //G4int ixy = x + (kModalitySize[1]-y-1)*kModalitySize[0];

	  G4int ixy = x + y*kModalitySize[0];
	  Index3D idx(x,y,z);
	  itr = kNestedModality.find(idx);
	  if(itr != kNestedModality.end()) {

	    modality[ixy] = kgMocrenIO->convertDensityToHU(itr->second);
	  } else {
	    modality[ixy] = -1024;
	  }

	}
      }
      kgMocrenIO->setModalityImage(modality);
    }

    //-- dose
    size_t nhits = kNestedHitsList.size();
    if(GFDEBUG) G4cout << "gMocren-file driver : # hits = " << nhits << G4endl;

    std::map<Index3D, G4double>::iterator hitsItr;
    std::map<G4String, std::map<Index3D, G4double> >::iterator hitsListItr = kNestedHitsList.begin();

    for(G4int n = 0; hitsListItr != kNestedHitsList.end(); hitsListItr++, n++) {

      kgMocrenIO->newDoseDist();
      kgMocrenIO->setDoseDistName(hitsListItr->first, n);
      kgMocrenIO->setDoseDistSize(kModalitySize, n);

      G4double minmax[2] = {DBL_MAX, -DBL_MAX};
      for(G4int z = 0 ; z < kModalitySize[2]; z++) {
	G4double * values = new G4double[nxy];
	for(G4int y = 0; y < kModalitySize[1]; y++) {
	  for(G4int x = 0; x < kModalitySize[0]; x++) {

	    G4int ixy = x + y*kModalitySize[0];
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
	kgMocrenIO->setDoseDist(values, n);
      }
      kgMocrenIO->setDoseDistMinMax(minmax, n);
      G4double lower = 0.;
      if(minmax[0] < 0)  lower = minmax[0];
      G4double scale = (minmax[1]-lower)/25000.;
      kgMocrenIO->setDoseDistScale(scale, n);
      G4String sunit("unit?"); //temporarily
      kgMocrenIO->setDoseDistUnit(sunit, n);
    }
      

    //-- draw axes
    if(false) {//true,false
      G4ThreeVector trans;
      G4RotationMatrix rot;
      trans = kVolumeTrans3D.getTranslation();
      rot = kVolumeTrans3D.getRotation().inverse();
      // x
      std::vector<G4float *> tracks;
      unsigned char colors[3];
      G4float * trk = new G4float[6];
      tracks.push_back(trk);

      G4ThreeVector orig(0.,0.,0), xa(2000.,0.,0.), ya(0.,2000.,0.), za(0.,0.,2000.);
      orig -= trans;
      orig.transform(rot);
      xa -= trans;
      xa.transform(rot);
      ya -= trans;
      ya.transform(rot);
      za -= trans;
      za.transform(rot);
      for(G4int i = 0; i < 3; i++) trk[i] = orig[i];
      for(G4int i = 0; i < 3; i++) trk[i+3] = xa[i];
      colors[0] = 255; colors[1] = 0; colors[2] = 0;
      kgMocrenIO->addTrack(tracks, colors);
      // y
      for(G4int i = 0; i < 3; i++) trk[i+3] = ya[i];
      colors[0] = 0; colors[1] = 255; colors[2] = 0;
      kgMocrenIO->addTrack(tracks, colors);
      // z
      for(G4int i = 0; i < 3; i++) trk[i+3] = za[i];
      colors[0] = 0; colors[1] = 0; colors[2] = 255;
      kgMocrenIO->addTrack(tracks, colors);
    }

    //-- detector
    ExtractDetector();


    if(GFDEBUG_DET) G4cout << ">>>>>>>>>>>>>>>>>>>>>>   (";
    std::vector<G4float> transformObjects;
    for(G4int i = 0; i < 3; i++) {
      // need to check!!
      transformObjects.push_back((kVolumeSize[i]/2. - kVoxelDimension[i]/2.));
      if(GFDEBUG_DET) G4cout << transformObjects[i] << ", ";
    }
    if(GFDEBUG_DET) G4cout << ")" << G4endl;


    kgMocrenIO->translateTracks(transformObjects);
    kgMocrenIO->translateDetector(transformObjects);

    // store
    kgMocrenIO->storeData(kGddFileName);
  } 

}


//----- 
void G4GMocrenFileSceneHandler::GFBeginModeling( void )
{
  G4VSceneHandler::BeginModeling();

  if( !GFIsInModeling() ) {

    if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "***** G4GMocrenFileSceneHandler::GFBeginModeling (called & started)" << G4endl;

    //----- Send saving command and heading comment
    BeginSavingGdd();

    kFlagInModeling = true ;

    // These models are entrusted to user commands /vis/scene/add/psHits or hits
    //GetScene()->AddEndOfEventModel(new G4PSHitsModel());
    //GetScene()->AddEndOfRunModel(new G4PSHitsModel());
    //scene->AddEndOfEventModel(new G4HitsModel());
    if(GFDEBUG_HIT) {
      G4Scene * scene = GetScene();
      std::vector<G4Scene::Model> vmodel = scene->GetEndOfEventModelList();
      std::vector<G4Scene::Model>::iterator itr = vmodel.begin();
      for(; itr != vmodel.end(); itr++) {
        G4cout << " IIIIII model name: " << itr->fpModel->GetGlobalTag() << G4endl;
      }
    }
  }
}


//========== AddPrimitive() functions ==========//

//----- Add polyline 
void G4GMocrenFileSceneHandler::AddPrimitive (const G4Polyline& polyline) 
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddPrimitive" << G4endl;

  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4GMocrenFileSceneHandler::AddPrimitive (const G4Polyline&)",
	 "gMocren1001", JustWarning,
	 "2D polylines not implemented.  Ignored.");
    }
    return;
  }

  //----- Initialize if necessary
  GFBeginModeling();

  static G4int numTrajectories = 0;
  if(numTrajectories >= MAX_NUM_TRAJECTORIES) return;

  // draw trajectories
  if(kbModelingTrajectory) {
    
    G4TrajectoriesModel * pTrModel = dynamic_cast<G4TrajectoriesModel*>(fpModel);
    if (!pTrModel) { 
      G4Exception 
	("G4VSceneHandler::AddCompound(const G4Polyline&)",
	 "gMocren0002", FatalException, "Not a G4TrajectoriesModel.");
    }

    G4ThreeVector trans;
    G4RotationMatrix rot;
    trans = kVolumeTrans3D.getTranslation();
    rot = kVolumeTrans3D.getRotation().inverse();

    if(GFDEBUG_TRK) G4cout << "   trajectory points : " << G4endl;
    std::vector<G4float *> trajectory;
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
      G4float * stepPts = new G4float[6];
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

    kgMocrenIO->addTrack(trajectory, trkcolor);

    numTrajectories++;
  }

} // G4GMocrenFileSceneHandler::AddPrimitive (polyline)


//----- Add text
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Text& text )
{
  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4GMocrenFileSceneHandler::AddPrimitive (const G4Text&)",
	 "gMocren1002", JustWarning,
	 "2D text not implemented.  Ignored.");
    }
    return;
  }

  // to avoid a warning in the compile process
  G4Text dummytext = text;

  //----- 
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddPrimitive( G4Text )" << G4endl;

  //----- Initialize IF NECESSARY
  GFBeginModeling();

} // G4GMocrenFileSceneHandler::AddPrimitive ( text )


//----- Add circle
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Circle& mark_circle )
{
  // to avoid a warning in the compile process
  G4Circle dummycircle = mark_circle;

  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4GMocrenFileSceneHandler::AddPrimitive (const G4Circle&)",
	 "gMocren1003", JustWarning,
	 "2D circles not implemented.  Ignored.");
    }
    return;
  }

  //----- 
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddPrimitive( G4Circle )" << G4endl;

  //----- Initialize IF NECESSARY
  GFBeginModeling();


} // G4GMocrenFileSceneHandler::AddPrimitive ( mark_circle )


//----- Add square
void G4GMocrenFileSceneHandler::AddPrimitive (const G4Square& mark_square )
{
  // to avoid a warning in the compile process
  G4Square dummysquare = mark_square;

  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4GMocrenFileSceneHandler::AddPrimitive (const G4Square&)",
	 "gMocren1004", JustWarning,
	 "2D squares not implemented.  Ignored.");
    }
    return;
  }

  //----- 
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddPrimitive( G4Square )" << G4endl;

  //----- Initialize if necessary
  GFBeginModeling();

} // G4GMocrenFileSceneHandler::AddPrimitive ( mark_square )


//----- Add polyhedron
void G4GMocrenFileSceneHandler::AddPrimitive ( const G4Polyhedron& polyhedron ) 
{
  //----- 
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddPrimitive( G4Polyhedron )" << G4endl;


  if (polyhedron.GetNoFacets() == 0) return;

  if (fProcessing2D) {
    static G4bool warned = false;
    if (!warned) {
      warned = true;
      G4Exception
	("G4GMocrenFileSceneHandler::AddPrimitive (const G4Polyhedron&)",
	 "gMocren1005", JustWarning,
	 "2D polyhedra not implemented.  Ignored.");
    }
    return;
  }

  //----- Initialize if necessary
  GFBeginModeling();

  //---------- (3) Facet block
  for (G4int f = polyhedron.GetNoFacets(); f; f--){
    G4bool notLastEdge = true;
    G4int index = -1; // initialization
    G4int edgeFlag = 1;
    //G4int preedgeFlag = 1;
    //G4int work[4], i = 0;
    G4int i = 0;
    do {
      //preedgeFlag = edgeFlag;
      notLastEdge = polyhedron.GetNextVertexIndex(index, edgeFlag);
      //work[i++] = index;
      i++;
    }while (notLastEdge);
    switch (i){
    case 3:
      //SendStrInt3(FR_FACET, work[0], work[1], work[2] );
      break;
    case 4:
      //SendStrInt4(FR_FACET, work[0], work[1], work[2], work[3] );
      break;
    default:
      if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
	G4cout <<
	  "ERROR G4GMocrenFileSceneHandler::AddPrimitive(G4Polyhedron)" << G4endl;
      G4PhysicalVolumeModel* pPVModel =
        dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
      if (pPVModel)   
	if(G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
	  G4cout << "Volume " << pPVModel->GetCurrentPV()->GetName() <<
	    ", Solid " << pPVModel->GetCurrentLV()->GetSolid()->GetName() <<
	    " (" << pPVModel->GetCurrentLV()->GetSolid()->GetEntityType();

      if(G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
	G4cout <<
	  "\nG4Polyhedron facet with " << i << " edges" << G4endl;	
    }
  }

} // G4GMocrenFileSceneHandler::AddPrimitive (polyhedron) 


//----- 
void G4GMocrenFileSceneHandler::GFEndModeling ()
{
  G4VSceneHandler::EndModeling();

  //----- 		
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** GFEndModeling (called)" << G4endl;

  if( GFIsInModeling() ) {

    if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations) {
      G4cout << "***** GFEndModeling (started) " ; 
      G4cout << "(/EndModeling, /DrawAll, /CloseDevice)" << G4endl;
    }

    //----- End saving data to g4.gdd
    EndSavingGdd() ;

    //------ Reset flag 
    kFlagInModeling = false ;

  }

}


//----- 
void G4GMocrenFileSceneHandler::BeginPrimitives (const G4Transform3D& objectTransformation)
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** BeginPrimitives " << G4endl;

  GFBeginModeling();

  G4VSceneHandler::BeginPrimitives (objectTransformation);

}


//----- 
void G4GMocrenFileSceneHandler::EndPrimitives ()
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** EndPrimitives " << G4endl;

  G4VSceneHandler::EndPrimitives ();
}


//========== AddSolid() functions ==========//

//----- Add box
void G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( box )" << G4endl;
  
  if(GFDEBUG_DET > 0)
    G4cout << "G4GMocrenFileSceneHandler::AddSolid(const G4Box&)  : "
	   << box.GetName() << G4endl;
  
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }
  
  //----- Initialize if necessary
  GFBeginModeling();
  
  
  //--
  if(GFDEBUG_DET > 1) {
    G4cout << "-------" << G4endl;
    G4cout << "    " << box.GetName() << G4endl;
    G4Polyhedron * poly = box.CreatePolyhedron();
    poly->Transform(fObjectTransformation);
    //G4int nv = poly->GetNoVertices();
    G4Point3D v1, v2;
    G4int next;
    //while(1) { // next flag isn't functional.
    for(G4int i = 0; i < 12; i++) { // # of edges is 12.
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
  
  
  // the volume name set by /vis/gMocren/setVolumeName
  G4String volName = kMessenger.getVolumeName();
  
  
  if(kFlagParameterization != 2) {
    G4ScoringManager * pScrMan = G4ScoringManager::GetScoringManager();
    if(pScrMan) {
      G4ScoringBox * pScBox = dynamic_cast<G4ScoringBox*>(pScrMan->FindMesh(volName));
      G4bool bMesh = false;
      if(pScBox != NULL) bMesh = true;
      if(bMesh) kFlagParameterization = 2;
      if(GFDEBUG_DET > 0) G4cout << "   G4ScoringManager::FindMesh() : "
        << volName << " - " << bMesh << G4endl;
    }
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
  
  //-- check the parameterised volume
  if(box.GetName() == volName) {
    
    kVolumeTrans3D = fObjectTransformation;
    // coordination system correction for gMocren
    G4ThreeVector raxis(1., 0., 0.), dummy(0.,0.,0.);
    G4RotationMatrix rot(raxis, pi*rad);
    G4Transform3D trot(rot, dummy);
    if(GFDEBUG_DET) {
      G4ThreeVector trans1 = kVolumeTrans3D.getTranslation();
      G4RotationMatrix rot1 = kVolumeTrans3D.getRotation().inverse();
      G4cout << "kVolumeTrans3D: " << trans1 << G4endl << rot1 << G4endl;
    }
    kVolumeTrans3D = kVolumeTrans3D*trot;
    if(GFDEBUG_DET) G4cout << " Parameterised volume : " << box.GetName() << G4endl;
    
    
    
    //
    G4VPhysicalVolume * pv[3] = {0,0,0};
    pv[0] = pPVModel->GetCurrentPV()->GetLogicalVolume()->GetDaughter(0);
    if(!pv[0]) {
      G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                  "gMocren0003", FatalException, "Unexpected volume.");
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
      default:
        G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                    "gMocren0004", FatalException, "Error.");
    }
    kNestedVolumeNames.push_back(pv[0]->GetName());
    if(GFDEBUG_DET)
      G4cout << "        daughter name :  " << pv[0]->GetName()
	     << "   # : " << nDaughters[0] << G4endl;
    
    //
    if(GFDEBUG_DET) {
      if(pv[0]->GetLogicalVolume()->GetNoDaughters()) {
        G4cout << "# of daughters : "
	       << pv[0]->GetLogicalVolume()->GetNoDaughters() << G4endl;
      } else {
        //G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
        //	    "gMocren0005", FatalException, "Error.");
      }
    }
    
    // check whether nested or regular parameterization
    if(GFDEBUG_DET) G4cout << "# of daughters : "
			   << pv[0]->GetLogicalVolume()->GetNoDaughters() << G4endl;
    if(pv[0]->GetLogicalVolume()->GetNoDaughters() == 0) {
      kFlagParameterization = 1;
      //G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
      //	    "gMocren0006", FatalException, "Error.");
    }
    
    if(kFlagParameterization == 0) {
      
      pv[1] = pv[0]->GetLogicalVolume()->GetDaughter(0);
      if(pv[1]) {
        pv[1]->GetReplicationData(axis, nReplicas, width, offset, consuming);
        nDaughters[1] = nReplicas;
        switch(axis) {
          case kXAxis: dirAxis[1] = 0; break;
          case kYAxis: dirAxis[1] = 1; break;
          case kZAxis: dirAxis[1] = 2; break;
          default:
            G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                        "gMocren0007", FatalException, "Error.");
        }
        kNestedVolumeNames.push_back(pv[1]->GetName());
        if(GFDEBUG_DET)
          G4cout << "        sub-daughter name :  " << pv[1]->GetName()
          << "   # : " << nDaughters[1]<< G4endl;
        
        //
        pv[2] = pv[1]->GetLogicalVolume()->GetDaughter(0);
        if(pv[2]) {
          nDaughters[2] = pv[2]->GetMultiplicity();
          kNestedVolumeNames.push_back(pv[2]->GetName());
          if(GFDEBUG_DET)
            G4cout << "        sub-sub-daughter name :  " << pv[2]->GetName()
            << "   # : " << nDaughters[2] << G4endl;
          
          if(nDaughters[2] > 1) {
            G4VNestedParameterisation * nestPara
            = dynamic_cast<G4VNestedParameterisation*>(pv[2]->GetParameterisation());
            if(nestPara == NULL)
              G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                          "gMocren0008", FatalException, "Non-nested parameterisation");
            
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
            else
              G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                          "gMocren0009", FatalException, "Unexpected nested parameterisation");
          }
        }
      }
      
      for(G4int i = 0; i < 3; i++) {
        kNestedVolumeDimension[i] = nDaughters[i];
        //kNestedVolumeDimension[i] = nDaughters[dirAxis[i]];
        kNestedVolumeDirAxis[i] = dirAxis[i];
      }
      //G4cout << "@@@@@@@@@ "
      //       << dirAxis[0] << ", " << dirAxis[1] << ", " << dirAxis[2] << G4endl;
      
      // get densities
      G4VNestedParameterisation * nestPara
      = dynamic_cast<G4VNestedParameterisation*>(pv[2]->GetParameterisation());
      if(nestPara != NULL) {
        G4double prexyz[3] = {0.,0.,0.}, xyz[3] = {0.,0.,0.};
        for(G4int n0 = 0; n0 < nDaughters[0]; n0++) {
          for(G4int n1 = 0; n1 < nDaughters[1]; n1++) {
            for(G4int n2 = 0; n2 < nDaughters[2]; n2++) {
              
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
                for(G4int i = 0; i < 3; i++) {
                  if(xyz[i] != prexyz[i])
                    G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                                "gMocren0010", FatalException, "Unsupported parameterisation");
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
              kNestedModality[i3d] = dens;
              if(GFDEBUG_DET)
                G4cout << " index: " << idx[0] << ", " << idx[1] << ", " << idx[2]
                << "  density: " << dens << G4endl;
              
              for(G4int i = 0; i < 3; i++) prexyz[i] = xyz[i];
            }
          }
        }
        
        kVolumeSize.set(box.GetXHalfLength()*2/mm,
                        box.GetYHalfLength()*2/mm,
                        box.GetZHalfLength()*2/mm);
        // mesh size
        if(!kbSetModalityVoxelSize) {
          G4float spacing[3] = {static_cast<G4float>(2*xyz[0]),
            static_cast<G4float>(2*xyz[1]),
            static_cast<G4float>(2*xyz[2])};
          kgMocrenIO->setVoxelSpacing(spacing);
          kVoxelDimension.set(spacing[0], spacing[1], spacing[2]);
          kbSetModalityVoxelSize = true;
        }
        
      } else {
        if(GFDEBUG_DET)
          G4cout << pv[2]->GetName() << G4endl;
        G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                    "gMocren0011", FatalException, "Non-nested parameterisation");
      }
      
      
      
      //-- debug
      if(GFDEBUG_DET > 1) {
        if(pPVModel->GetCurrentPV()->IsParameterised()) {
          G4VPVParameterisation * para = pPVModel->GetCurrentPV()->GetParameterisation();
          G4cout << " Is nested parameterisation? : " << para->IsNested() << G4endl;
          
          
          G4int npvp = pPVModel->GetDrawnPVPath().size();
          G4cout << "     physical volume node id : "
          << "size: " << npvp << ", PV name: ";
          for(G4int i = 0; i < npvp; i++) {
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
          
          
          pPVModel->GetCurrentPV()->GetReplicationData(axis, nReplicas, width, offset, consuming);
          G4cout << "     # replicas : " << nReplicas << G4endl;
          G4double pareDims[3] = {0.,0.,0.};
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
      
      
    } else if(kFlagParameterization == 1) { // G4PhantomParameterisation based geom. construnction
      
      // get the dimension of the parameterized patient geometry
      G4PhantomParameterisation * phantomPara
      = dynamic_cast<G4PhantomParameterisation*>(pv[0]->GetParameterisation());
      if(phantomPara == NULL) {
        G4Exception("G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )",
                    "gMocren0012", FatalException, "no G4PhantomParameterisation");
      } else {
        ;
      }
      
      kNestedVolumeDimension[0] = phantomPara->GetNoVoxelX();
      kNestedVolumeDimension[1] = phantomPara->GetNoVoxelY();
      kNestedVolumeDimension[2] = phantomPara->GetNoVoxelZ();
      kNestedVolumeDirAxis[0] = 0;
      kNestedVolumeDirAxis[1] = 1;
      kNestedVolumeDirAxis[2] = 2;
      
      // get densities of the parameterized patient geometry
      G4int nX = kNestedVolumeDimension[0];
      G4int nXY = kNestedVolumeDimension[0]*kNestedVolumeDimension[1];
      
      for(G4int n0 = 0; n0 < kNestedVolumeDimension[0]; n0++) {
        for(G4int n1 = 0; n1 < kNestedVolumeDimension[1]; n1++) {
          for(G4int n2 = 0; n2 < kNestedVolumeDimension[2]; n2++) {
            
            G4int repNo = n0 + n1*nX + n2*nXY;
            G4Material * mat = phantomPara->ComputeMaterial(repNo, pv[0]);
            G4double dens = mat->GetDensity()/(g/cm3);
            
            
            G4int idx[3];
            idx[kNestedVolumeDirAxis[0]] = n0;
            idx[kNestedVolumeDirAxis[1]] = n1;
            idx[kNestedVolumeDirAxis[2]] = n2;
            Index3D i3d(idx[0],idx[1],idx[2]);
            kNestedModality[i3d] = dens;
            
            if(GFDEBUG_DET)
              G4cout << " index: " << idx[0] << ", " << idx[1] << ", " << idx[2]
              << "  density: " << dens << G4endl;
            
          }
        }
      }
      
      kVolumeSize.set(box.GetXHalfLength()*2/mm,
                      box.GetYHalfLength()*2/mm,
                      box.GetZHalfLength()*2/mm);
      
      // mesh size
      if(!kbSetModalityVoxelSize) {
        G4float spacing[3] = {static_cast<G4float>(2*phantomPara->GetVoxelHalfX()),
          static_cast<G4float>(2*phantomPara->GetVoxelHalfY()),
          static_cast<G4float>(2*phantomPara->GetVoxelHalfZ())};
        kgMocrenIO->setVoxelSpacing(spacing);
        kVoxelDimension.set(spacing[0], spacing[1], spacing[2]);
        kbSetModalityVoxelSize = true;
      }
    }
    
  } // if(box.GetName() == volName)
  
  
  // processing geometry construction based on the interactive PS
  if(!kFlagProcessedInteractiveScorer) {
    
    
    // get the dimension of the geometry defined in G4VScoringMesh
    G4ScoringManager * pScrMan = G4ScoringManager::GetScoringManager();
    //if(!pScrMan) return;
    if(pScrMan) {
    G4ScoringBox * scoringBox
    = dynamic_cast<G4ScoringBox*>(pScrMan->FindMesh(volName));
    //if(scoringBox == NULL) return;
    if(scoringBox) {

    
    
    G4int nVoxels[3];
    scoringBox->GetNumberOfSegments(nVoxels);
    // this order depends on the G4ScoringBox
    kNestedVolumeDimension[0] = nVoxels[2];
    kNestedVolumeDimension[1] = nVoxels[1];
    kNestedVolumeDimension[2] = nVoxels[0];
    kNestedVolumeDirAxis[0] = 2;
    kNestedVolumeDirAxis[1] = 1;
    kNestedVolumeDirAxis[2] = 0;
    
    // get densities of the parameterized patient geometry
    for(G4int n0 = 0; n0 < kNestedVolumeDimension[0]; n0++) {
      for(G4int n1 = 0; n1 < kNestedVolumeDimension[1]; n1++) {
        for(G4int n2 = 0; n2 < kNestedVolumeDimension[2]; n2++) {
          
          G4double dens = 0.*(g/cm3);
          
          G4int idx[3];
          idx[kNestedVolumeDirAxis[0]] = n0;
          idx[kNestedVolumeDirAxis[1]] = n1;
          idx[kNestedVolumeDirAxis[2]] = n2;
          Index3D i3d(idx[0],idx[1],idx[2]);
          kNestedModality[i3d] = dens;
          
        }
      }
    }
    
    G4ThreeVector boxSize = scoringBox->GetSize();
    if(GFDEBUG_DET > 1) {
      G4cout << "Interactive Scorer : size - "
	     << boxSize.x()/cm << " x "
	     << boxSize.y()/cm << " x "
	     << boxSize.z()/cm << " [cm3]" << G4endl;
      G4cout << "Interactive Scorer : # voxels - "
	     << nVoxels[0] << " x "
	     << nVoxels[1] << " x "
	     << nVoxels[2] << G4endl;
    }
    kVolumeSize.set(boxSize.x()*2,
                    boxSize.y()*2,
                    boxSize.z()*2);
    
    // mesh size
    if(!kbSetModalityVoxelSize) {
      G4float spacing[3] = {static_cast<G4float>(boxSize.x()*2/nVoxels[0]),
        static_cast<G4float>(boxSize.y()*2/nVoxels[1]),
        static_cast<G4float>(boxSize.z()*2/nVoxels[2])};
      
      kgMocrenIO->setVoxelSpacing(spacing);
      kVoxelDimension.set(spacing[0], spacing[1], spacing[2]);
      kbSetModalityVoxelSize = true;
      
    }
    
    
    kVolumeTrans3D = fObjectTransformation;
    
    // translation for the scoring mesh
    G4ThreeVector sbth = scoringBox->GetTranslation();
    G4Translate3D sbtranslate(sbth);
    kVolumeTrans3D = kVolumeTrans3D*sbtranslate;
    
    // rotation matrix for the scoring mesh 
    G4RotationMatrix sbrm;
    sbrm = scoringBox->GetRotationMatrix();
    if(!sbrm.isIdentity()) {
      G4ThreeVector sbdummy(0.,0.,0.); 
      G4Transform3D sbrotate(sbrm.inverse(), sbdummy);
      kVolumeTrans3D = kVolumeTrans3D*sbrotate;
    }
    
    
    // coordination system correction for gMocren
    G4ThreeVector raxisY(0., 1., 0.), dummyY(0.,0.,0.); 
    G4RotationMatrix rotY(raxisY, pi*rad);
    G4Transform3D trotY(rotY, dummyY);
    G4ThreeVector raxisZ(0., 0., 1.), dummyZ(0.,0.,0.); 
    G4RotationMatrix rotZ(raxisZ, pi*rad);
    G4Transform3D trotZ(rotZ, dummyZ);
    
    kVolumeTrans3D = kVolumeTrans3D*trotY*trotZ;
    
    }
    }
    //
    kFlagProcessedInteractiveScorer = true;
  }
  
  
  static G4VPhysicalVolume * volPV = NULL;
  if(pPVModel->GetCurrentPV()->GetName() == volName) {
    volPV = pPVModel->GetCurrentPV();
  }
  
  //-- add detectors
  G4bool bAddDet = true;
  if(!kMessenger.getDrawVolumeGrid()) {

    if(kFlagParameterization == 0) { // nested parameterisation
      /*
      G4String volDSolidName;
      if(volPV) {
        G4int nDaughter = volPV->GetLogicalVolume()->GetNoDaughters();
        G4VPhysicalVolume * volDPV = NULL;
        if(nDaughter > 0) volDPV = volPV->GetLogicalVolume()->GetDaughter(0);
        if(volDPV) {
          nDaughter = volDPV->GetLogicalVolume()->GetNoDaughters();
          if(nDaughter > 0)
            volDSolidName = volDPV->GetLogicalVolume()->GetDaughter(0)
                           ->GetLogicalVolume()->GetSolid()->GetName();
        }
      }
      */
      
      //std::cout << "Parameterization volume: " << volName << " - "
      //          << box.GetName() << std::endl;
      
      if(volName == box.GetName()) {
        bAddDet = false;
      }
      
      std::vector<G4String>::iterator itr = kNestedVolumeNames.begin();
      for(; itr != kNestedVolumeNames.end(); itr++) {
        if(*itr == box.GetName())  {
          bAddDet = false;
          break;
        }
      }
    } else if(kFlagParameterization == 1) { // phantom paramemterisation
      
      G4String volDSolidName;
      if(volPV) {
        volDSolidName = volPV->GetLogicalVolume()->GetDaughter(0)
                        ->GetLogicalVolume()->GetSolid()->GetName();
      }

      //std::cout << "Phantom Parameterization volume: " << volDSolidName
      //          << " - " << box.GetName() << std::endl;
      
      if(volDSolidName == box.GetName()) {
        bAddDet = false;
      }
      
    } else if(kFlagParameterization == 2) { // interactive primitive scorer
      //std::cout << "Regular Parameterization volume: " << box.GetName() << std::endl;
    }
    
  }
  if(bAddDet) AddDetector(box);
  

} // void G4GMocrenFileSceneHandler::AddSolid( const G4Box& box )


//----- Add tubes
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Tubs& tubes )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( tubes )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(tubes);


  // for a debug
  if(GFDEBUG_DET > 0) {
    G4cout << "-------" << G4endl;
    G4cout << "    " << tubes.GetName() << G4endl;
    G4Polyhedron * poly = tubes.CreatePolyhedron();
    G4int nv = poly->GetNoVertices();
    for(G4int i = 0; i < nv; i++) {
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

} // void G4GMocrenFileSceneHandler::AddSolid( const G4Tubs& )



//----- Add cons
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Cons& cons )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( cons )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(cons);

}// G4GMocrenFileSceneHandler::AddSolid( cons )


//----- Add trd
void G4GMocrenFileSceneHandler::AddSolid ( const G4Trd& trd )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( trd )" << G4endl;


  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(trd);

} // G4GMocrenFileSceneHandler::AddSolid ( trd )


//----- Add sphere
void G4GMocrenFileSceneHandler::AddSolid ( const G4Sphere& sphere )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( sphere )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(sphere);

} // G4GMocrenFileSceneHandler::AddSolid ( sphere )


//----- Add para
void G4GMocrenFileSceneHandler::AddSolid (const G4Para& para)
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( para )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(para);

} // G4GMocrenFileSceneHandler::AddSolid ( para )


//----- Add trap
void G4GMocrenFileSceneHandler::AddSolid (const G4Trap& trap)
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( trap )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(trap);

} // G4GMocrenFileSceneHandler::AddSolid (const G4Trap& trap)


//----- Add torus
void 
G4GMocrenFileSceneHandler::AddSolid( const G4Torus& torus )
{
  if(GFDEBUG || G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "***** AddSolid ( torus )" << G4endl;

  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

  //
  AddDetector(torus);

} // void G4GMocrenFileSceneHandler::AddSolid( const G4Torus& )



//----- Add a shape which is not treated above
void G4GMocrenFileSceneHandler::AddSolid ( const G4VSolid& solid  )
{
  //----- skip drawing invisible primitive
  if( !IsVisible() ) { return ; }

  //----- Initialize if necessary
  GFBeginModeling();

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

  kbModelingTrajectory = true;

  G4VSceneHandler::AddCompound(traj);

  if(GFDEBUG_TRK) {
    G4cout << " ::AddCompound(const G4VTrajectory&) >>>>>>>>> " << G4endl;
    G4TrajectoriesModel * pTrModel = dynamic_cast<G4TrajectoriesModel*>(fpModel);
    if (!pTrModel) { 
      G4Exception 
	("G4VSceneHandler::AddCompound(const G4VTrajectory&)",
	 "gMocren0013", FatalException, "Not a G4TrajectoriesModel.");
    } else {
      traj.DrawTrajectory();

      const G4VTrajectory * trj = pTrModel->GetCurrentTrajectory();
      G4cout << "------ track" << G4endl;
      G4cout << "    name:     " << trj->GetParticleName() << G4endl;
      G4cout << "    id:       " << trj->GetTrackID() << G4endl;
      G4cout << "    charge:   " << trj->GetCharge() << G4endl;
      G4cout << "    momentum: " << trj->GetInitialMomentum() << G4endl;
      
      G4int nPnt = trj->GetPointEntries();
      G4cout << "    point:    ";
      for(G4int i = 0; i < nPnt; i++) {
	G4cout << trj->GetPoint(i)->GetPosition() << ", ";
      }
      G4cout << G4endl;
    }
    G4cout << G4endl;
  }

  kbModelingTrajectory = false;
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

  std::vector<G4String> hitNames = kMessenger.getHitNames();
  if(GFDEBUG_HIT) {
    std::vector<G4String>::iterator itr = hitNames.begin();
    for(; itr != hitNames.end(); itr++) 
      G4cout << "  hit name : " << *itr << G4endl;
  }
  
  std::vector<G4AttValue> * attval = hit.CreateAttValues();
  if(!attval) {G4cout << "0 empty " << G4endl;}
  else {

    G4bool bid[3] = {false, false, false};
    Index3D id;

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

    G4int nhitname = (G4int)hitNames.size();

    if(bid[0] && bid[1] && bid[2]) {

      if(GFDEBUG_HIT)
	G4cout << " Hit : index(" << id.x << ", " << id.y << ", "
	       << id.z << ")" << G4endl;

      // Get attributes
      for(itr = attval->begin(); itr != attval->end(); itr++) {
	for(G4int i = 0; i < nhitname; i++) {
	  if(itr->GetName() == hitNames[i]) {

	    std::string stmp = itr->GetValue();
	    std::istringstream sval(stmp.c_str());
	    G4double value;
	    G4String unit;
	    sval >> value >> unit;

	    std::map<G4String, std::map<Index3D, G4double> >::iterator kNestedHitsListItr;
	    kNestedHitsListItr = kNestedHitsList.find(hitNames[i]);
	    if(kNestedHitsListItr != kNestedHitsList.end()) {
	      //fTempNestedHits = &kNestedHitsListItr->second;
	      //(*fTempNestedHits)[id] = value;
	      kNestedHitsListItr->second[id] = value;
	    } else {
	      std::map<Index3D, G4double> hits;
	      hits.insert(std::map<Index3D, G4double>::value_type(id, value));
	      kNestedHitsList[hitNames[i]] = hits;
	    }

	    
	    if(GFDEBUG_HIT)
	      G4cout << "     : " << hitNames[i] << " -> " << value
		     << " [" << unit << "]" << G4endl;
	  }
	}
      }
    } else {
      G4Exception("G4GMocrenFileSceneHandler::AddCompound(const G4VHit &)",
		  "gMocren0014", FatalException, "Error");
    }

    delete attval;
  }

}

void G4GMocrenFileSceneHandler::AddCompound( const G4VDigi & digi) {
  if(GFDEBUG_DIGI) G4cout << " ::AddCompound(const G4VDigi&) >>>>>>>>> " << G4endl;
  G4VSceneHandler::AddCompound(digi);
}

void G4GMocrenFileSceneHandler::AddCompound(const G4THitsMap<G4double> & hits) {
  if(GFDEBUG_HIT)
    G4cout << " ::AddCompound(const std::map<G4int, G4double*> &) >>>>>>>>> " << G4endl;


  std::vector<G4String> hitScorerNames = kMessenger.getHitScorerNames();
  G4int nhitname = (G4int)hitScorerNames.size();
  G4String scorername = static_cast<G4VHitsCollection>(hits).GetName();

  //-- --//
  /*
  std::vector<G4String> hitScorerNames = kMessenger.getHitScorerNames();
  if(GFDEBUG_HIT) {
    std::vector<G4String>::iterator itr = hitScorerNames.begin();
    for(; itr != hitScorerNames.end(); itr++) 
      G4cout << "  PS name : " << *itr << G4endl;
  }
  */
  
  {  // Scope bracket to avoid compiler messages about shadowing (JA).
  //for(G4int i = 0; i < nhitname; i++) {       // this selection trusts
    //if(scorername == hitScorerNames[i]) {   // thea command /vis/scene/add/psHits hit_name.

      G4int idx[3];
      std::map<G4int, G4double*> * map = hits.GetMap();
      std::map<G4int, G4double*>::const_iterator itr = map->begin();
      for(; itr != map->end(); itr++) {
	GetNestedVolumeIndex(itr->first, idx);
	Index3D id(idx[0], idx[1], idx[2]);
	
	std::map<G4String, std::map<Index3D, G4double> >::iterator nestedHitsListItr;
	nestedHitsListItr = kNestedHitsList.find(scorername);
	if(nestedHitsListItr != kNestedHitsList.end()) {
	  nestedHitsListItr->second[id] = *(itr->second);
	} else {
	  std::map<Index3D, G4double> hit;
	  hit.insert(std::map<Index3D, G4double>::value_type(id, *(itr->second)));
	  kNestedHitsList[scorername] = hit;
	}
      }
 
      //break;
    //}
  //}
  }

  if(GFDEBUG_HIT) {
    G4String meshname = static_cast<G4VHitsCollection>(hits).GetSDname();
    G4cout << "       >>>>> " << meshname << " : " << scorername  << G4endl;

    for(G4int i = 0; i < nhitname; i++)
      if(scorername == hitScorerNames[i]) 
	G4cout << "       !!!! Hit scorer !!!! " << scorername << G4endl;

    G4cout << " dimension: "
	   << kNestedVolumeDimension[0] << " x "
	   << kNestedVolumeDimension[1] << " x "
	   << kNestedVolumeDimension[2] << G4endl;

    G4int id[3];
    std::map<G4int, G4double*> * map = hits.GetMap();
    std::map<G4int, G4double*>::const_iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      GetNestedVolumeIndex(itr->first, id);
      G4cout << "[" << itr->first << "] "
	     << "("<< id[0] << "," << id[1] << "," << id[2] << ")"
	     << *(itr->second) << ", ";
    }
    G4cout << G4endl;
  }
}

void G4GMocrenFileSceneHandler::AddCompound(const G4THitsMap<G4StatDouble> & hits) {
  if(GFDEBUG_HIT)
    G4cout << " ::AddCompound(const std::map<G4int, G4StatDouble*> &) >>>>>>>>> " << G4endl;


  std::vector<G4String> hitScorerNames = kMessenger.getHitScorerNames();
  G4int nhitname = (G4int)hitScorerNames.size();
  G4String scorername = static_cast<G4VHitsCollection>(hits).GetName();

  //-- --//
  /*
  std::vector<G4String> hitScorerNames = kMessenger.getHitScorerNames();
  if(GFDEBUG_HIT) {
    std::vector<G4String>::iterator itr = hitScorerNames.begin();
    for(; itr != hitScorerNames.end(); itr++) 
      G4cout << "  PS name : " << *itr << G4endl;
  }
  */
  
  {  // Scope bracket to avoid compiler messages about shadowing (JA).
  //for(G4int i = 0; i < nhitname; i++) {       // this selection trusts
    //if(scorername == hitScorerNames[i]) {   // thea command /vis/scene/add/psHits hit_name.

      G4int idx[3];
      std::map<G4int, G4StatDouble*> * map = hits.GetMap();
      std::map<G4int, G4StatDouble*>::const_iterator itr = map->begin();
      for(; itr != map->end(); itr++) {
	GetNestedVolumeIndex(itr->first, idx);
	Index3D id(idx[0], idx[1], idx[2]);
	
	std::map<G4String, std::map<Index3D, G4double> >::iterator nestedHitsListItr;
	nestedHitsListItr = kNestedHitsList.find(scorername);
	if(nestedHitsListItr != kNestedHitsList.end()) {
	  nestedHitsListItr->second[id] = itr->second->sum_wx();
	} else {
	  std::map<Index3D, G4double> hit;
	  hit.insert(std::map<Index3D, G4double>::value_type(id, itr->second->sum_wx()));
	  kNestedHitsList[scorername] = hit;
	}
      }
 
      //break;
    //}
  //}
  }

  if(GFDEBUG_HIT) {
    G4String meshname = static_cast<G4VHitsCollection>(hits).GetSDname();
    G4cout << "       >>>>> " << meshname << " : " << scorername  << G4endl;

    for(G4int i = 0; i < nhitname; i++)
      if(scorername == hitScorerNames[i]) 
	G4cout << "       !!!! Hit scorer !!!! " << scorername << G4endl;

    G4cout << " dimension: "
	   << kNestedVolumeDimension[0] << " x "
	   << kNestedVolumeDimension[1] << " x "
	   << kNestedVolumeDimension[2] << G4endl;

    G4int id[3];
    std::map<G4int, G4StatDouble*> * map = hits.GetMap();
    std::map<G4int, G4StatDouble*>::const_iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      GetNestedVolumeIndex(itr->first, id);
      G4cout << "[" << itr->first << "] "
	     << "("<< id[0] << "," << id[1] << "," << id[2] << ")"
	     << itr->second->sum_wx() << ", ";
    }
    G4cout << G4endl;
  }
}

//----- 
G4bool G4GMocrenFileSceneHandler::IsVisible()
{
  //----- 
  G4bool  visibility  = true ;

  const G4VisAttributes* pVisAttribs =
    fpViewer->GetApplicableVisAttributes( fpVisAttribs );

  if(pVisAttribs) {
    visibility = pVisAttribs->IsVisible();
  } 

  return visibility ;

} // G4GMocrenFileSceneHandler::IsVisible()


//----- 
void G4GMocrenFileSceneHandler::ClearTransientStore() 
{
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
  std::vector<G4float *> dedges;
  G4Polyhedron * poly = solid.CreatePolyhedron();
  detector.polyhedron = poly;
  detector.transform3D = fObjectTransformation;

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
  for(G4int i = 0; i < 3; i++) detector.color[i] = uccolor[i];
  //
  kDetectors.push_back(detector);

  if(GFDEBUG_DET > 1) {
    G4cout << "0     color:   (" << (G4int)uccolor[0] << ", "
	   << (G4int)uccolor[1] << ", " << (G4int)uccolor[2] << ")"
	   << G4endl;
  }

}

//----- 
void G4GMocrenFileSceneHandler::ExtractDetector() {

  std::vector<Detector>::iterator itr = kDetectors.begin();

  for(; itr != kDetectors.end(); itr++) {

    // detector name
    G4String detname = itr->name;
    if(GFDEBUG_DET > 1)
      G4cout << "Detector name : " << detname << G4endl;

    // edge points of the detector
    std::vector<G4float *> dedges;
    G4Polyhedron * poly = itr->polyhedron;
    poly->Transform(itr->transform3D);
    G4Transform3D invVolTrans = kVolumeTrans3D.inverse();
    poly->Transform(invVolTrans);

    G4Point3D v1, v2;
    G4bool bnext = true;
    G4int next;
    G4int nedges = 0;
    //
    while(bnext) {
      if(!(poly->GetNextEdge(v1, v2, next))) bnext =false;
      G4float * edge = new G4float[6];
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
    kgMocrenIO->addDetector(detname, dedges, uccolor);
    for(G4int i = 0; i < nedges; i++) { // # of edges is 12.
      delete [] dedges[i];
    }
    dedges.clear(); 

    if(GFDEBUG_DET > 1) {
      G4cout << "    color:   (" << (G4int)uccolor[0] << ", "
	     << (G4int)uccolor[1] << ", " << (G4int)uccolor[2] << ")"
	     << G4endl;
    }
  }
}

void G4GMocrenFileSceneHandler::GetNestedVolumeIndex(G4int _idx, G4int _idx3d[3]) {
  if(kNestedVolumeDimension[0] == 0 ||
     kNestedVolumeDimension[1] == 0 ||
     kNestedVolumeDimension[2] == 0) {
    for(G4int i = 0; i < 3; i++) _idx3d[i] = 0;
    return;
  }


  if(kFlagParameterization == 0) {

    G4int plane = kNestedVolumeDimension[2]*kNestedVolumeDimension[1];
    G4int line = kNestedVolumeDimension[2];
 
  /*
  G4int idx3d[3];
  idx3d[0] = _idx/plane;
  idx3d[1] = (_idx%plane)/line;
  idx3d[2] = (_idx%plane)%line;
  _idx3d[0] = idx3d[kNestedVolumeDirAxis[0]];
  _idx3d[1] = idx3d[kNestedVolumeDirAxis[1]];
  _idx3d[2] = idx3d[kNestedVolumeDirAxis[2]];
  */

    _idx3d[kNestedVolumeDirAxis[0]] = _idx/plane;
    _idx3d[kNestedVolumeDirAxis[1]] = (_idx%plane)/line;
    _idx3d[kNestedVolumeDirAxis[2]] = (_idx%plane)%line;



  /*

  G4cout << "G4GMocrenFileSceneHandler::GetNestedVolumeIndex : " << G4endl;
  G4cout << "(depi, depj, depk) : "
	 << kNestedVolumeDirAxis[0] << ", "
	 << kNestedVolumeDirAxis[1] << ", "
	 << kNestedVolumeDirAxis[2] << G4endl;
  G4cout << "(ni, nj, nk) :"
	 << kNestedVolumeDimension[0] << ", " 
	 << kNestedVolumeDimension[1] << ", "
	 << kNestedVolumeDimension[2] << " - " << G4endl;

  G4cout << " _idx = " << _idx << "  :  plane = "
	 << plane << " ,   line = " << line << G4endl;
  G4cout << "(idx,idy,idz) + " << _idx3d[0] << ", "
	 << _idx3d[1] << ", " << _idx3d[2] << " + " << G4endl;

  */



  } else {
    
    G4int plane = kNestedVolumeDimension[0]*kNestedVolumeDimension[1];
    G4int line = kNestedVolumeDimension[0];
    _idx3d[kNestedVolumeDirAxis[2]] = _idx/plane;
    _idx3d[kNestedVolumeDirAxis[1]] = (_idx%plane)/line;
    _idx3d[kNestedVolumeDirAxis[0]] = (_idx%plane)%line;

  }

}


//-- --//
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

//-- --//
G4GMocrenFileSceneHandler::Index3D::Index3D()
  : x(0), y(0), z(0) {
  ;
}

G4GMocrenFileSceneHandler::Index3D::Index3D(const Index3D & _index3D) 
  : x(_index3D.x), y(_index3D.y), z(_index3D.z) {
  //: x(_index3D.X()),
  //y(_index3D.Y()),
  //z(_index3D.Z()) {
  //  : x(static_cast<Index3D>(_index3D).x),
  //    y(static_cast<Index3D>(_index3D).y),
  //    z(static_cast<Index3D>(_index3D).z) {
  ;
}

G4GMocrenFileSceneHandler::Index3D::Index3D(G4int _x, G4int _y, G4int _z) 
  : x(_x), y(_y), z(_z) {
  ;
}
G4bool G4GMocrenFileSceneHandler::Index3D::operator < (const Index3D & _right) const {
  if(z < static_cast<Index3D>(_right).z) {
     return true;
  } else if(z == _right.z) {
    if(y < static_cast<Index3D>(_right).y) return true;
    else if(y == _right.y) 
      if(x < static_cast<Index3D>(_right).x) return true;
  } 
  return false;
}
G4bool G4GMocrenFileSceneHandler::Index3D::operator == (const Index3D & _right) const {
  if(z == _right.z && y == _right.y && x == _right.x) return true;
  return false;
}
