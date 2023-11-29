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
//
//
// File I/O manager class for writing or reading calcuated dose
// distribution and some event information
//
// Created:  Mar. 31, 2009  Akinori Kimura : release for the gMocrenFile driver
//
//                               Akinori Kimura
//                               gMocren home page:
//                               http://geant4.kek.jp/gMocren/
//
//
#include "G4GMocrenIO.hh"
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

#include "globals.hh"
#include "G4VisManager.hh"

#if defined(_WIN32)
#define LITTLE_ENDIAN 1234
#define BYTE_ORDER LITTLE_ENDIAN
#endif

const int DOSERANGE = 25000;

//----- GMocrenDataPrimitive class in the GMocrenDataIO class-----//
template <typename T> 
GMocrenDataPrimitive<T>::GMocrenDataPrimitive () {
  clear();
}
template <typename T> 
GMocrenDataPrimitive<T>::~GMocrenDataPrimitive () {
  /*
    std::vector<short *>::iterator itr = image.begin();
    for(; itr != image.end(); itr++) {
    delete [] *itr;
    }
  */
}

template <typename T> GMocrenDataPrimitive<T> & 
GMocrenDataPrimitive<T>::operator = (const GMocrenDataPrimitive<T> & _right) {
  if (this == &_right) return *this;
  for(int i = 0; i < 3; i++) {
    kSize[i] = _right.kSize[i];
    kCenter[i] = _right.kCenter[i];
  }
  kScale = _right.kScale;
  for(int i = 0; i < 2; i++) kMinmax[i] = _right.kMinmax[i];
  int num = kSize[0]*kSize[1];
  kImage.clear();
  for(int z = 0; z < kSize[2]; z++) {
    T * img = new T[num];
    for(int i = 0; i < num; i++) img[i] =_right.kImage[z][i];
    kImage.push_back(img);
  }
  return *this;
}

template <typename T> GMocrenDataPrimitive<T> & 
GMocrenDataPrimitive<T>::operator + (const GMocrenDataPrimitive<T> & _right) {

  GMocrenDataPrimitive<T> rprim;
  bool stat = true;
  for(int i = 0; i < 3; i++) {
    if(kSize[i] != _right.kSize[i]) stat = false;
    if(kCenter[i] != _right.kCenter[i]) stat = false;
  }
  if(!stat) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Warning: operator + "
	     << "         Cannot do the operator +"
	     << G4endl;
    return *this;
  }

  rprim.setSize(kSize);
  rprim.setCenterPosition(kCenter);
  
  T mms[2] = {9e100,-9e100};
  //if(mms[0] > _right.minmax[0]) mms[0] = _right.minmax[0];
  //if(mms[1] < _right.minmax[1]) mms[1] = _right.minmax[1];

  int num = kSize[0]*kSize[1];
  for(int z = 0; z < kSize[2]; z++) {
    T * img = new T[num];
    for(int xy = 0; xy < num; xy++) {
      img[xy] = kImage[z][xy] + _right.kImage[z][xy];
      if(mms[0] > img[xy]) mms[0] = img[xy];
      if(mms[1] < img[xy]) mms[1] = img[xy];
    }
    rprim.addImage(img);
  }
  rprim.setMinMax(mms);

  T scl = mms[1]/DOSERANGE;
  rprim.setScale(scl);

  return rprim;
}

template <typename T> GMocrenDataPrimitive<T> & 
GMocrenDataPrimitive<T>::operator += (const GMocrenDataPrimitive<T> & _right) {

  bool stat = true;
  for(int i = 0; i < 3; i++) {
    if(kSize[i] != _right.kSize[i]) stat = false;
    if(kCenter[i] != _right.kCenter[i]) stat = false;
  }
  if(!stat) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Warning: operator += " << G4endl
	     << "         Cannot do the operator +="
	     << G4endl;
    return *this;
  }

  if(kMinmax[0] > _right.kMinmax[0]) kMinmax[0] = _right.kMinmax[0];
  if(kMinmax[1] < _right.kMinmax[1]) kMinmax[1] = _right.kMinmax[1];

  int num = kSize[0]*kSize[1];
  for(int z = 0; z < kSize[2]; z++) {
    for(int xy = 0; xy < num; xy++) {
      kImage[z][xy] += _right.kImage[z][xy];
      if(kMinmax[0] > kImage[z][xy]) kMinmax[0] = kImage[z][xy];
      if(kMinmax[1] < kImage[z][xy]) kMinmax[1] = kImage[z][xy];
    }
  }

  kScale = kMinmax[1]/DOSERANGE;

  return *this;
}

template <typename T> 
void GMocrenDataPrimitive<T>::clear() {
  for(int i = 0; i < 3; i++) {
    kSize[i] = 0;
    kCenter[i] = 0.;
  }
  kScale = 1.;
  kMinmax[0] = (T)32109;
  kMinmax[1] = (T)-32109;

  clearImage();
}
template <typename T> 
void GMocrenDataPrimitive<T>::clearImage() {
  typename std::vector<T *>::iterator itr;
  for(itr = kImage.begin(); itr != kImage.end(); itr++) {
    delete [] *itr;
  }
  kImage.clear();
}
template <typename T> 
void GMocrenDataPrimitive<T>::setSize(int _size[3]) {
  for(int i = 0; i < 3; i++) kSize[i] = _size[i];
}
template <typename T> 
void GMocrenDataPrimitive<T>::getSize(int _size[3]) {
  for(int i = 0; i < 3; i++) _size[i] = kSize[i];
}
template <typename T> 
void GMocrenDataPrimitive<T>::setScale(double & _scale) {
  kScale = _scale;
}
template <typename T> 
double GMocrenDataPrimitive<T>::getScale() {
  return kScale;
}
template <typename T> 
void GMocrenDataPrimitive<T>::setMinMax(T _minmax[2]) {
  for(int i = 0; i < 2; i++) kMinmax[i] = _minmax[i];
}
template <typename T> 
void GMocrenDataPrimitive<T>::getMinMax(T _minmax[2]) {
  for(int i = 0; i < 2; i++) _minmax[i] = kMinmax[i];

}
template <typename T> 
void GMocrenDataPrimitive<T>::setImage(std::vector<T *> & _image) {
  kImage = _image;
}
template <typename T> 
void GMocrenDataPrimitive<T>::addImage(T * _image) {
  kImage.push_back(_image);
}
template <typename T> 
std::vector<T *> & GMocrenDataPrimitive<T>::getImage() {
  return kImage;
}
template <typename T> 
T * GMocrenDataPrimitive<T>::getImage(int _z) {
  if(_z >= (int)kImage.size())  return 0;
  return kImage[_z];
}
template <typename T> 
void GMocrenDataPrimitive<T>::setCenterPosition(float _center[3]) {
  for(int i = 0; i < 3; i++) kCenter[i] = _center[i];
}
template <typename T> 
void GMocrenDataPrimitive<T>::getCenterPosition(float _center[3]) {
  for(int i = 0; i < 3; i++) _center[i] = kCenter[i];
}
template <typename T> 
void GMocrenDataPrimitive<T>::setName(std::string & _name) {
  kDataName = _name;
}
template <typename T> 
std::string GMocrenDataPrimitive<T>::getName() {
  return kDataName;
}





GMocrenTrack::GMocrenTrack() {
    kTrack.clear();
    for(int i = 0; i < 3; i++) kColor[i] = 0;
}

void GMocrenTrack::addStep(float _startx, float _starty, float _startz,
			   float _endx, float _endy, float _endz) {
  struct Step step;
  step.startPoint[0] = _startx;
  step.startPoint[1] = _starty;
  step.startPoint[2] = _startz;
  step.endPoint[0] = _endx;
  step.endPoint[1] = _endy;
  step.endPoint[2] = _endz;
  kTrack.push_back(step);
}
void GMocrenTrack::getStep(float & _startx, float & _starty, float & _startz,
			   float & _endx, float & _endy, float & _endz,
			   int _num) {
  if(_num >= (int)kTrack.size()) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "GMocrenTrack::getStep(...) Error: "
	     << "invalid step # : " << _num << G4endl;
    return;
  }

  _startx = kTrack[_num].startPoint[0];
  _starty = kTrack[_num].startPoint[1];
  _startz = kTrack[_num].startPoint[2];
  _endx = kTrack[_num].endPoint[0];
  _endy = kTrack[_num].endPoint[1];
  _endz = kTrack[_num].endPoint[2];
}
void GMocrenTrack::translate(std::vector<float> & _translate) {
  std::vector<struct Step>::iterator itr = kTrack.begin();
  for(; itr != kTrack.end(); itr++) {
    for(int i = 0; i < 3; i++ ) {
      itr->startPoint[i] += _translate[i];
      itr->endPoint[i] += _translate[i];
    }
  } 
}









GMocrenDetector::GMocrenDetector() {
    kDetector.clear();
    for(int i = 0; i < 3; i++) kColor[i] = 0;
}

void GMocrenDetector::addEdge(float _startx, float _starty, float _startz,
			      float _endx, float _endy, float _endz) {
  struct Edge edge;
  edge.startPoint[0] = _startx;
  edge.startPoint[1] = _starty;
  edge.startPoint[2] = _startz;
  edge.endPoint[0] = _endx;
  edge.endPoint[1] = _endy;
  edge.endPoint[2] = _endz;
  kDetector.push_back(edge);
}
void GMocrenDetector::getEdge(float & _startx, float & _starty, float & _startz,
			   float & _endx, float & _endy, float & _endz,
			   int _num) {
  if(_num >= (int)kDetector.size()) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "GMocrenDetector::getEdge(...) Error: "
	     << "invalid edge # : " << _num << G4endl;
    return;
  }

  _startx = kDetector[_num].startPoint[0];
  _starty = kDetector[_num].startPoint[1];
  _startz = kDetector[_num].startPoint[2];
  _endx = kDetector[_num].endPoint[0];
  _endy = kDetector[_num].endPoint[1];
  _endz = kDetector[_num].endPoint[2];
}
void GMocrenDetector::translate(std::vector<float> & _translate) {
  std::vector<struct Edge>::iterator itr = kDetector.begin();
  for(; itr != kDetector.end(); itr++) {
    for(int i = 0; i < 3; i++) {
      itr->startPoint[i] += _translate[i];
      itr->endPoint[i] += _translate[i];
    } 
  }
}









// file information
std::string G4GMocrenIO::kId;
std::string G4GMocrenIO::kVersion = "2.0.0";
int G4GMocrenIO::kNumberOfEvents = 0;
char G4GMocrenIO::kLittleEndianInput = true;

#if BYTE_ORDER == LITTLE_ENDIAN
char G4GMocrenIO::kLittleEndianOutput = true;
#else
char G4GMocrenIO::kLittleEndianOutput = false; // Big endian
#endif
std::string G4GMocrenIO::kComment;
//
std::string G4GMocrenIO::kFileName = "dose.gdd";

//
unsigned int G4GMocrenIO::kPointerToModalityData = 0;
std::vector<unsigned int> G4GMocrenIO::kPointerToDoseDistData;
unsigned int G4GMocrenIO::kPointerToROIData = 0;
unsigned int G4GMocrenIO::kPointerToTrackData = 0;
unsigned int G4GMocrenIO::kPointerToDetectorData = 0;

// modality
float G4GMocrenIO::kVoxelSpacing[3] = {0., 0., 0.};
class GMocrenDataPrimitive<short>  G4GMocrenIO::kModality;
std::vector<float> G4GMocrenIO::kModalityImageDensityMap;
std::string G4GMocrenIO::kModalityUnit = "g/cm3       "; // 12 Bytes

// dose
std::vector<class GMocrenDataPrimitive<double> > G4GMocrenIO::kDose;
std::string G4GMocrenIO::kDoseUnit = "keV         "; // 12 Bytes

// ROI
std::vector<class GMocrenDataPrimitive<short> > G4GMocrenIO::kRoi;

// track
std::vector<float *> G4GMocrenIO::kSteps;
std::vector<unsigned char *> G4GMocrenIO::kStepColors;
std::vector<class GMocrenTrack> G4GMocrenIO::kTracks;

// detector
std::vector<class GMocrenDetector> G4GMocrenIO::kDetectors;

// verbose
int G4GMocrenIO::kVerbose = 0;

const int IDLENGTH  = 21;
const int VERLENGTH = 6;

// constructor
G4GMocrenIO::G4GMocrenIO()
  : kTracksWillBeStored(true) {
  ;
}

// destructor
G4GMocrenIO::~G4GMocrenIO() {
  ;
}

// initialize
void G4GMocrenIO::initialize() {

  kId.clear();
  kVersion = "2.0.0";
  kNumberOfEvents = 0;
  kLittleEndianInput = true;
#if BYTE_ORDER == LITTLE_ENDIAN
  kLittleEndianOutput = true;
#else // Big endian
  kLittleEndianOutput = false;
#endif
  kComment.clear();
  kFileName = "dose.gdd";
  kPointerToModalityData = 0;
  kPointerToDoseDistData.clear();
  kPointerToROIData = 0;
  kPointerToTrackData = 0;
  // modality
  for(int i = 0; i < 3; i++) kVoxelSpacing[i] = 0.;
  kModality.clear();
  kModalityImageDensityMap.clear();
  kModalityUnit = "g/cm3       "; // 12 Bytes
  // dose
  kDose.clear();
  kDoseUnit = "keV         "; // 12 Bytes
  // ROI
  kRoi.clear();
  // track
  std::vector<float *>::iterator itr;
  for(itr = kSteps.begin(); itr != kSteps.end(); itr++) delete [] *itr;
  kSteps.clear();
  std::vector<unsigned char *>::iterator citr;
  for(citr = kStepColors.begin(); citr != kStepColors.end(); citr++)
    delete [] *citr;
  kStepColors.clear();
  kTracksWillBeStored = true;

  // verbose
  kVerbose = 0;
}

bool G4GMocrenIO::storeData() {
  return storeData4();
}
//
bool G4GMocrenIO::storeData(char * _filename) {
  return storeData4(_filename);
}

bool G4GMocrenIO::storeData4() {

  bool DEBUG = false;//

  if(DEBUG || kVerbose > 0)
    G4cout << ">>>>>>>  store data (ver.4) <<<<<<<" << G4endl;
  if(DEBUG || kVerbose > 0)
    G4cout << "         " << kFileName << G4endl;

  // output file open
  std::ofstream ofile(kFileName.c_str(),
		      std::ios_base::out|std::ios_base::binary);
  if(DEBUG || kVerbose > 0)
    G4cout << "         file open status: " << ofile.rdbuf() << G4endl;
  
  // file identifier
  ofile.write("gMocren ", 8);

  // file version
  unsigned char ver = 0x04;
  ofile.write((char *)&ver, 1);

  // endian
  //ofile.write((char *)&kLittleEndianOutput, sizeof(char));
  char littleEndian = 0x01;
  ofile.write((char *)&littleEndian, sizeof(char));
  if(DEBUG || kVerbose > 0) {
    //G4cout << "Endian: " << (int)kLittleEndianOutput << G4endl;
    G4cout << "Endian: " << (int)littleEndian << G4endl;
  }

  // for inverting the byte order
  float ftmp[6];
  int itmp[6];
  short stmp[6];

  // comment length (fixed size)
  int commentLength = 1024;
  if(kLittleEndianOutput) {
    ofile.write((char *)&commentLength, 4);
  } else {
    invertByteOrder((char *)&commentLength, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  // comment 
  char cmt[1025];
  std::strncpy(cmt, kComment.c_str(), 1024);
  cmt[1024] = '\0';
  ofile.write(cmt, 1024);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Data comment : "
	      << kComment << G4endl;
  }

  // voxel spacings for all images
  if(kLittleEndianOutput) {
    ofile.write((char *)kVoxelSpacing, 12);
  } else {
    for(int j = 0; j < 3; j++)
      invertByteOrder((char *)&kVoxelSpacing[j], ftmp[j]);
    ofile.write((char *)ftmp, 12);
  }
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }

  calcPointers4();
  if(!kTracksWillBeStored) kPointerToTrackData = 0;

  // offset from file starting point to the modality image data
  if(kLittleEndianOutput) {
    ofile.write((char *)&kPointerToModalityData, 4);
  } else {
    invertByteOrder((char *)&kPointerToModalityData, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  // # of dose distributions
  //int nDoseDist = (int)pointerToDoseDistData.size();
  int nDoseDist = getNumDoseDist();
  if(kLittleEndianOutput) {
    ofile.write((char *)&nDoseDist, 4);
  } else {
    invertByteOrder((char *)&nDoseDist, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  // offset from file starting point to the dose image data
  if(kLittleEndianOutput) {
    for(int i = 0; i < nDoseDist; i++) {
      ofile.write((char *)&kPointerToDoseDistData[i], 4);
    }
  } else {
    for(int i = 0; i < nDoseDist; i++) {
      invertByteOrder((char *)&kPointerToDoseDistData[i], itmp[0]);
      ofile.write((char *)itmp, 4);
    }
  }

  // offset from file starting point to the ROI image data
  if(kLittleEndianOutput) {
    ofile.write((char *)&kPointerToROIData, 4);
  } else {
    invertByteOrder((char *)&kPointerToROIData, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  // offset from file starting point to the track data
  if(kLittleEndianOutput) {
    ofile.write((char *)&kPointerToTrackData, 4);
  } else {
    invertByteOrder((char *)&kPointerToTrackData, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  // offset from file starting point to the detector data
  if(kLittleEndianOutput) {
    ofile.write((char *)&kPointerToDetectorData, 4);
  } else {
    invertByteOrder((char *)&kPointerToDetectorData, itmp[0]);
    ofile.write((char *)itmp, 4);
  }

  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", ";
    for(int i = 0; i < nDoseDist; i++) {
      G4cout << kPointerToDoseDistData[i] << ", ";
    }
    G4cout << kPointerToROIData << ", "
	      << kPointerToTrackData << ", "
	      << kPointerToDetectorData
	      << G4endl;
  }

  //----- modality image -----//

  int size[3];
  float scale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];
  // modality image size
  kModality.getSize(size);

  if(kLittleEndianOutput) {
    ofile.write((char *)size, 3*sizeof(int));
  } else {
    for(int j = 0; j < 3; j++)
      invertByteOrder((char *)&size[j], itmp[j]);
    ofile.write((char *)itmp, 12);
  }

  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }

  // modality image max. & min.
  kModality.getMinMax(minmax);
  if(kLittleEndianOutput) {
    ofile.write((char *)minmax, 4);
  } else {
    for(int j = 0; j < 2; j++)
      invertByteOrder((char *)&minmax[j], stmp[j]);
    ofile.write((char *)stmp, 4);
  }

  // modality image unit
  char munit[13] = "g/cm3\0";
  ofile.write((char *)munit, 12);

  // modality image scale
  scale = (float)kModality.getScale();
  if(kLittleEndianOutput) {
    ofile.write((char *)&scale, 4);
  } else {
    invertByteOrder((char *)&scale, ftmp[0]);
    ofile.write((char *)ftmp, 4);
  }
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image min., max., scale : "
	      << minmax[0] << ", "
	      << minmax[1] << ", "
	      << scale << G4endl;
  }

  // modality image
  int psize = size[0]*size[1];
  if(DEBUG || kVerbose > 0) G4cout << "Modality image : ";
  for(int i = 0; i < size[2]; i++) {
    short * image = kModality.getImage(i);
    if(kLittleEndianOutput) {
      ofile.write((char *)image, psize*sizeof(short));
    } else {
      for(int j = 0; j < psize; j++) {
	invertByteOrder((char *)&image[j], stmp[0]);
	ofile.write((char *)stmp, 2);
      }
    }

    if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << image[(size_t)(psize*0.55)] << ", ";
  }
  if(DEBUG || kVerbose > 0) G4cout << G4endl;

  // modality desity map for CT value
  size_t msize = minmax[1] - minmax[0]+1;
  if(DEBUG || kVerbose > 0) 
    G4cout << "modality image : " << minmax[0] << ", " << minmax[1] << G4endl;
  float * pdmap = new float[msize];
  for(int i = 0; i < (int)msize; i++) pdmap[i] =kModalityImageDensityMap[i]; 

  if(kLittleEndianOutput) {
    ofile.write((char *)pdmap, msize*sizeof(float));
  } else {
    for(int j = 0; j < (int)msize; j++) {
      invertByteOrder((char *)&pdmap[j], ftmp[0]);
      ofile.write((char *)ftmp, 4);
    }
  }

  if(DEBUG || kVerbose > 0) {
    G4cout << "density map : " << std::ends;
    for(int i = 0; i < (int)msize; i+=50)
      G4cout <<kModalityImageDensityMap[i] << ", ";
    G4cout << G4endl;
  }
  delete [] pdmap;


  //----- dose distribution image -----//

  if(!isDoseEmpty()) {

    calcDoseDistScale();

    for(int ndose = 0; ndose < nDoseDist; ndose++) {
      // dose distrbution image size
      kDose[ndose].getSize(size);
      if(kLittleEndianOutput) {
	ofile.write((char *)size, 3*sizeof(int));
      } else {
	for(int j = 0; j < 3; j++)
	  invertByteOrder((char *)&size[j], itmp[j]);
	ofile.write((char *)itmp, 12);
      }
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose << "] image size : ("
		  << size[0] << ", "
		  << size[1] << ", "
		  << size[2] << ")"
		  << G4endl;
      }

      // dose distribution max. & min.
      getShortDoseDistMinMax(minmax, ndose);
      if(kLittleEndianOutput) {
	ofile.write((char *)minmax, 2*2); // sizeof(shorft)*2
      } else {
	for(int j = 0; j < 2; j++)
	  invertByteOrder((char *)&minmax[j], stmp[j]);
	ofile.write((char *)stmp, 4);
      }

      // dose distribution unit
      char cdunit[13];
      std::strncpy(cdunit, kDoseUnit.c_str(), 12);
      cdunit[12] = '\0';
      ofile.write(cdunit, 12);
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. unit : " << kDoseUnit << G4endl;
      }

      // dose distribution scaling 
      double dscale;
      dscale = getDoseDistScale(ndose);
      scale = float(dscale);
      if(kLittleEndianOutput) {
	ofile.write((char *)&scale, 4);
      } else {
	invertByteOrder((char *)&scale, ftmp[0]);
	ofile.write((char *)ftmp, 4);
      }
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose
		  << "] image min., max., scale : "
		  << minmax[0] << ", "
		  << minmax[1] << ", "
		  << scale << G4endl;
      }

      // dose distribution image
      int dsize = size[0]*size[1];
      short * dimage = new short[dsize];
      for(int z = 0; z < size[2]; z++) {
	getShortDoseDist(dimage, z, ndose);
	if(kLittleEndianOutput) {
	  ofile.write((char *)dimage, dsize*2); //sizeof(short)
	} else {
	  for(int j = 0; j < dsize; j++) {
	    invertByteOrder((char *)&dimage[j], stmp[0]);
	    ofile.write((char *)stmp, 2);
	  }
	}

	if(DEBUG || kVerbose > 0) {
	  for(int j = 0; j < dsize; j++) {
	    if(dimage[j] < 0)
	      G4cout << "[" << j << "," << z << "]"
			<< dimage[j] << ", ";
	  }
	}
      }
      if(DEBUG || kVerbose > 0) G4cout << G4endl;
      delete [] dimage;

      // relative location of the dose distribution image for 
      // the modality image
      getDoseDistCenterPosition(fCenter, ndose);
      for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
      if(kLittleEndianOutput) {
	ofile.write((char *)iCenter, 3*4); // 3*sizeof(int)
      } else {
	for(int j = 0; j < 3; j++)
	  invertByteOrder((char *)&iCenter[j], itmp[j]);
	ofile.write((char *)itmp, 12);
      }
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose
		  << "]image relative location : ("
		  << iCenter[0] << ", "
		  << iCenter[1] << ", "
		  << iCenter[2] << ")" << G4endl;
      }

      // dose distribution name
      std::string name = getDoseDistName(ndose);
      if(name.size() == 0) name = "dose";
      name.resize(80);
      ofile.write((char *)name.c_str(), 80);
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. name : " << name << G4endl;
      }

    }
  }

  //----- ROI image -----//
  if(!isROIEmpty()) {
    // ROI image size
    kRoi[0].getSize(size);
    if(kLittleEndianOutput) {
      ofile.write((char *)size, 3*sizeof(int));
    } else {
      for(int j = 0; j < 3; j++)
	invertByteOrder((char *)&size[j], itmp[j]);
      ofile.write((char *)itmp, 12);
    }
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    kRoi[0].getMinMax(minmax);
    if(kLittleEndianOutput) {
      ofile.write((char *)minmax, sizeof(short)*2);
    } else {
      for(int j = 0; j < 2; j++)
	invertByteOrder((char *)&minmax[j], stmp[j]);
      ofile.write((char *)stmp, 4);
    }

    // ROI distribution scaling 
    scale = (float)kRoi[0].getScale();
    if(kLittleEndianOutput) {
      ofile.write((char *)&scale, sizeof(float));
    } else {
      invertByteOrder((char *)&scale, ftmp[0]);
      ofile.write((char *)ftmp, 4);
    }
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    for(int i = 0; i < size[2]; i++) {
      short * rimage = kRoi[0].getImage(i);
      if(kLittleEndianOutput) {
	ofile.write((char *)rimage, rsize*sizeof(short));
      } else {
	for(int j = 0; j < rsize; j++) {
	  invertByteOrder((char *)&rimage[j], stmp[0]);
	  ofile.write((char *)stmp, 2);
	}
      }

    }

    // ROI relative location
    kRoi[0].getCenterPosition(fCenter);
    for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
    if(kLittleEndianOutput) {
      ofile.write((char *)iCenter, 3*sizeof(int));
    } else {
      for(int j = 0; j < 3; j++)
	invertByteOrder((char *)&iCenter[j], itmp[j]);
      ofile.write((char *)itmp, 12);
    }
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< iCenter[0] << ", "
		<< iCenter[1] << ", "
		<< iCenter[2] << ")" << G4endl;
    }
  }

  //----- track information -----//
  // number of track 
  if(kPointerToTrackData > 0) {

    int ntrk = (int)kTracks.size();
    if(kLittleEndianOutput) {
      ofile.write((char *)&ntrk, sizeof(int));
    } else {
      invertByteOrder((char *)&ntrk, itmp[0]);
      ofile.write((char *)itmp, 4);
    }
    if(DEBUG || kVerbose > 0) {
      G4cout << "# of tracks : "
		<< ntrk << G4endl;
    }

    for(int nt = 0; nt < ntrk; nt++) {

      // # of steps in a track
      int nsteps = kTracks[nt].getNumberOfSteps();
      if(kLittleEndianOutput) {
	ofile.write((char *)&nsteps, sizeof(int));
      } else {
	invertByteOrder((char *)&nsteps, itmp[0]);
	ofile.write((char *)itmp, 4);
      }
      if(DEBUG || kVerbose > 0) {
	G4cout << "# of steps : " << nsteps << G4endl;
      }

      // track color
      unsigned char tcolor[3];
      kTracks[nt].getColor(tcolor);
      ofile.write((char *)tcolor, 3);

      // steps
      float stepPoints[6];
      for(int isteps = 0; isteps < nsteps; isteps++) {
	kTracks[nt].getStep(stepPoints[0], stepPoints[1], stepPoints[2],
			    stepPoints[3], stepPoints[4], stepPoints[5],
			    isteps);

	if(kLittleEndianOutput) {
	  ofile.write((char *)stepPoints, sizeof(float)*6);
	} else {
	  for(int j = 0; j < 6; j++)
	    invertByteOrder((char *)&stepPoints[j], ftmp[j]);
	  ofile.write((char *)ftmp, 24);
	}
      }
    }
  }

  //----- detector information -----//
  // number of detectors
  if(kPointerToDetectorData > 0) {
    int ndet = (int)kDetectors.size();
    if(kLittleEndianOutput) {
      ofile.write((char *)&ndet, sizeof(int));
    } else {
      invertByteOrder((char *)&ndet, itmp[0]);
      ofile.write((char *)itmp, 4);
    }
    if(DEBUG || kVerbose > 0) {
      G4cout << "# of detectors : "
		<< ndet << G4endl;
    }

    for(int nd = 0; nd < ndet; nd++) {

      // # of edges of a detector
      int nedges = kDetectors[nd].getNumberOfEdges();
      if(kLittleEndianOutput) {
	ofile.write((char *)&nedges, sizeof(int));
      } else {
	invertByteOrder((char *)&nedges, itmp[0]);
	ofile.write((char *)itmp, 4);
      }
      if(DEBUG || kVerbose > 0) {
	G4cout << "# of edges in a detector : " << nedges << G4endl;
      }

      // edges
      float edgePoints[6];
      for(int ne = 0; ne < nedges; ne++) {
	kDetectors[nd].getEdge(edgePoints[0], edgePoints[1], edgePoints[2],
			       edgePoints[3], edgePoints[4], edgePoints[5],
			       ne);

	if(kLittleEndianOutput) {
	  ofile.write((char *)edgePoints, sizeof(float)*6);
	} else {
	  for(int j = 0; j < 6; j++)
	    invertByteOrder((char *)&edgePoints[j], ftmp[j]);
	  ofile.write((char *)ftmp, 24);
	}

	if(DEBUG || kVerbose > 0) {
	  if(ne < 1) {
	    G4cout << " edge : (" << edgePoints[0] << ", "
		      << edgePoints[1] << ", "
		      << edgePoints[2] << ") - ("
		      << edgePoints[3] << ", "
		      << edgePoints[4] << ", "
		      << edgePoints[5] << ")" << G4endl;
	  }
	}
      }

      // detector color
      unsigned char dcolor[3];
      kDetectors[nd].getColor(dcolor);
      ofile.write((char *)dcolor, 3);
      if(DEBUG || kVerbose > 0) {
	G4cout << " rgb : (" << (int)dcolor[0] << ", "
		  << (int)dcolor[1] << ", "
		  << (int)dcolor[2] << ")" << G4endl;
      }

      // detector name
      std::string dname = kDetectors[nd].getName();
      dname.resize(80);
      ofile.write((char *)dname.c_str(), 80);
      if(DEBUG || kVerbose > 0) {
	G4cout << " detector name : " << dname << G4endl;
      
      }
    }
  }

  // file end mark
  ofile.write("END", 3);

  ofile.close();
  if(DEBUG || kVerbose > 0)
    G4cout << ">>>> closed gdd file: " << kFileName << G4endl;

  return true;
}
bool G4GMocrenIO::storeData3() {

  if(kVerbose > 0) G4cout << ">>>>>>>  store data (ver.3) <<<<<<<" << G4endl;
  if(kVerbose > 0) G4cout << "         " << kFileName << G4endl;

  bool DEBUG = false;//

  // output file open
  std::ofstream ofile(kFileName.c_str(),
		      std::ios_base::out|std::ios_base::binary);

  // file identifier
  ofile.write("gMocren ", 8);

  // file version
  unsigned char ver = 0x03;
  ofile.write((char *)&ver, 1);

  // endian
  ofile.write((char *)&kLittleEndianOutput, sizeof(char));

  // comment length (fixed size)
  int commentLength = 1024;
  ofile.write((char *)&commentLength, 4);

  // comment 
  char cmt[1025];
  std::strncpy(cmt, kComment.c_str(), 1024);
  ofile.write((char *)cmt, 1024);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Data comment : "
	      << kComment << G4endl;
  }

  // voxel spacings for all images
  ofile.write((char *)kVoxelSpacing, 12);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }

  calcPointers3();

  // offset from file starting point to the modality image data
  ofile.write((char *)&kPointerToModalityData, 4);

  // # of dose distributions
  //int nDoseDist = (int)pointerToDoseDistData.size();
  int nDoseDist = getNumDoseDist();
  ofile.write((char *)&nDoseDist, 4);

  // offset from file starting point to the dose image data
  for(int i = 0; i < nDoseDist; i++) {
    ofile.write((char *)&kPointerToDoseDistData[i], 4);
  }

  // offset from file starting point to the ROI image data
  ofile.write((char *)&kPointerToROIData, 4);

  // offset from file starting point to the track data
  ofile.write((char *)&kPointerToTrackData, 4);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", ";
    for(int i = 0; i < nDoseDist; i++) {
      G4cout << kPointerToDoseDistData[i] << ", ";
    }
    G4cout << kPointerToROIData << ", "
	      << kPointerToTrackData << G4endl;
  }

  //----- modality image -----//

  int size[3];
  float scale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];
  // modality image size
  kModality.getSize(size);
  ofile.write((char *)size, 3*sizeof(int));
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }

  // modality image max. & min.
  kModality.getMinMax(minmax);
  ofile.write((char *)minmax, 4);

  // modality image unit
  char munit[13] = "g/cm3       ";
  ofile.write((char *)munit, 12);

  // modality image scale
  scale = (float)kModality.getScale();
  ofile.write((char *)&scale, 4);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image min., max., scale : "
	      << minmax[0] << ", "
	      << minmax[1] << ", "
	      << scale << G4endl;
  }

  // modality image
  int psize = size[0]*size[1];
  if(DEBUG || kVerbose > 0) G4cout << "Modality image : ";
  for(int i = 0; i < size[2]; i++) {
    short * image = kModality.getImage(i);
    ofile.write((char *)image, psize*sizeof(short));

    if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << image[(size_t)(psize*0.55)] << ", ";
  }
  if(DEBUG || kVerbose > 0) G4cout << G4endl;

  // modality desity map for CT value
  size_t msize = minmax[1] - minmax[0]+1;
  float * pdmap = new float[msize];
  for(int i = 0; i < (int)msize; i++) pdmap[i] =kModalityImageDensityMap[i]; 
  ofile.write((char *)pdmap, msize*sizeof(float));
  if(DEBUG || kVerbose > 0) {
    G4cout << "density map : " << std::ends;
    for(int i = 0; i < (int)msize; i+=50)
      G4cout <<kModalityImageDensityMap[i] << ", ";
    G4cout << G4endl;
  }
  delete [] pdmap;


  //----- dose distribution image -----//

  if(!isDoseEmpty()) {

    calcDoseDistScale();

    for(int ndose = 0; ndose < nDoseDist; ndose++) {
      // dose distrbution image size
      kDose[ndose].getSize(size);
      ofile.write((char *)size, 3*sizeof(int));
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose << "] image size : ("
		  << size[0] << ", "
		  << size[1] << ", "
		  << size[2] << ")"
		  << G4endl;
      }

      // dose distribution max. & min.
      getShortDoseDistMinMax(minmax, ndose);
      ofile.write((char *)minmax, 2*2); // sizeof(shorft)*2

      // dose distribution unit
      ofile.write((char *)kDoseUnit.c_str(), 12);
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. unit : " << kDoseUnit << G4endl;
      }

      // dose distribution scaling 
      double dscale;
      dscale = getDoseDistScale(ndose);
      scale = float(dscale);
      ofile.write((char *)&scale, 4);
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose
		  << "] image min., max., scale : "
		  << minmax[0] << ", "
		  << minmax[1] << ", "
		  << scale << G4endl;
      }

      // dose distribution image
      int dsize = size[0]*size[1];
      short * dimage = new short[dsize];
      for(int z = 0; z < size[2]; z++) {
	getShortDoseDist(dimage, z, ndose);
	ofile.write((char *)dimage, dsize*2); //sizeof(short)

	if(DEBUG || kVerbose > 0) {
	  for(int j = 0; j < dsize; j++) {
	    if(dimage[j] < 0)
	      G4cout << "[" << j << "," << z << "]"
			<< dimage[j] << ", ";
	  }
	}
      }
      if(DEBUG || kVerbose > 0) G4cout << G4endl;
      delete [] dimage;

      // relative location of the dose distribution image for 
      // the modality image
      getDoseDistCenterPosition(fCenter, ndose);
      for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
      ofile.write((char *)iCenter, 3*4); // 3*sizeof(int)
      if(DEBUG || kVerbose > 0) {
	G4cout << "Dose dist. [" << ndose
		  << "]image relative location : ("
		  << iCenter[0] << ", "
		  << iCenter[1] << ", "
		  << iCenter[2] << ")" << G4endl;
      }
    }
  }

  //----- ROI image -----//
  if(!isROIEmpty()) {
    // ROI image size
    kRoi[0].getSize(size);
    ofile.write((char *)size, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    kRoi[0].getMinMax(minmax);
    ofile.write((char *)minmax, sizeof(short)*2);

    // ROI distribution scaling 
    scale = (float)kRoi[0].getScale();
    ofile.write((char *)&scale, sizeof(float));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    for(int i = 0; i < size[2]; i++) {
      short * rimage = kRoi[0].getImage(i);
      ofile.write((char *)rimage, rsize*sizeof(short));

    }

    // ROI relative location
    kRoi[0].getCenterPosition(fCenter);
    for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
    ofile.write((char *)iCenter, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< iCenter[0] << ", "
		<< iCenter[1] << ", "
		<< iCenter[2] << ")" << G4endl;
    }
  }

  //----- track information -----//
  // number of track 
  int ntrk = (int)kSteps.size();
  ofile.write((char *)&ntrk, sizeof(int));
  if(DEBUG || kVerbose > 0) {
    G4cout << "# of tracks : "
	      << ntrk << G4endl;
  }
  // track position
  for(int i = 0; i < ntrk; i++) {
    float * tp = kSteps[i];
    ofile.write((char *)tp, sizeof(float)*6);
  }
  // track color
  int ntcolor = int(kStepColors.size());
  if(ntrk != ntcolor) 
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "# of track color information must be the same as # of tracks." 
	     << G4endl;
  unsigned char white[3] = {255,255,255}; // default color
  for(int i = 0; i < ntrk; i++) {
    if(i < ntcolor) {
      unsigned char * tcolor = kStepColors[i];
      ofile.write((char *)tcolor, 3);
    } else {
      ofile.write((char *)white, 3);
    }
  }
  
  // file end mark
  ofile.write("END", 3);

  ofile.close();

  return true;
}
//
bool G4GMocrenIO::storeData4(char * _filename) {
  kFileName = _filename;
  return storeData4();
}

// version 2
bool G4GMocrenIO::storeData2() {

  if(kVerbose > 0) G4cout << ">>>>>>>  store data (ver.2) <<<<<<<" << G4endl;
  if(kVerbose > 0) G4cout << "         " << kFileName << G4endl;

  bool DEBUG = false;//

  // output file open
  std::ofstream ofile(kFileName.c_str(),
		      std::ios_base::out|std::ios_base::binary);

  // file identifier
  ofile.write("GRAPE    ", 8);

  // file version
  unsigned char ver = 0x02;
  ofile.write((char *)&ver, 1);
  // file id for old file format support
  ofile.write(kId.c_str(), IDLENGTH);
  // file version for old file format support
  ofile.write(kVersion.c_str(), VERLENGTH);
  // endian
  ofile.write((char *)&kLittleEndianOutput, sizeof(char));

  /*
  // event number
  ofile.write((char *)&numberOfEvents, sizeof(int));
  float imageSpacing[3]; 
  imageSpacing[0] = modalityImageVoxelSpacing[0];
  imageSpacing[1] = modalityImageVoxelSpacing[1];
  imageSpacing[2] = modalityImageVoxelSpacing[2];
  ofile.write((char *)imageSpacing, 12);
  */


  // voxel spacings for all images
  ofile.write((char *)kVoxelSpacing, 12);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }

  calcPointers2();
  // offset from file starting point to the modality image data
  ofile.write((char *)&kPointerToModalityData, 4);

  // offset from file starting point to the dose image data
  ofile.write((char *)&kPointerToDoseDistData[0], 4);

  // offset from file starting point to the ROI image data
  ofile.write((char *)&kPointerToROIData, 4);

  // offset from file starting point to the track data
  ofile.write((char *)&kPointerToTrackData, 4);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", "
	      << kPointerToDoseDistData[0] << ", "
	      << kPointerToROIData << ", "
	      << kPointerToTrackData << G4endl;
  }

  //----- modality image -----//

  int size[3];
  float scale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];
  // modality image size
  kModality.getSize(size);
  ofile.write((char *)size, 3*sizeof(int));
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }

  // modality image max. & min.
  kModality.getMinMax(minmax);
  ofile.write((char *)minmax, 4);

  // modality image unit
  //char munit[13] = "g/cm3       ";
  //ofile.write((char *)&munit, 12);
  
  // modality image scale
  scale = (float)kModality.getScale();
  ofile.write((char *)&scale, 4);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image min., max., scale : "
	      << minmax[0] << ", "
	      << minmax[1] << ", "
	      << scale << G4endl;
  }

  // modality image
  int psize = size[0]*size[1];
  if(DEBUG || kVerbose > 0) G4cout << "Modality image : ";
  for(int i = 0; i < size[2]; i++) {
    short * image =kModality.getImage(i);
    ofile.write((char *)image, psize*sizeof(short));

    if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << image[(size_t)(psize*0.55)] << ", ";
  }
  if(DEBUG || kVerbose > 0) G4cout << G4endl;

  // modality desity map for CT value
  size_t msize = minmax[1] - minmax[0]+1;
  float * pdmap = new float[msize];
  for(int i = 0; i < (int)msize; i++) pdmap[i] =kModalityImageDensityMap[i]; 
  ofile.write((char *)pdmap, msize*sizeof(float));
  if(DEBUG || kVerbose > 0) {
    G4cout << "density map : " << std::ends;
    for(int i = 0; i < (int)msize; i+=50)
      G4cout <<kModalityImageDensityMap[i] << ", ";
    G4cout << G4endl;
  }
  delete [] pdmap;


  //----- dose distribution image -----//

  if(!isDoseEmpty()) {
    calcDoseDistScale();

    // dose distrbution image size
    kDose[0].getSize(size);
    ofile.write((char *)size, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // dose distribution max. & min.
    getShortDoseDistMinMax(minmax);
    ofile.write((char *)minmax, sizeof(short)*2);

    // dose distribution scaling 
    scale = (float)kDose[0].getScale();
    ofile.write((char *)&scale, sizeof(float));
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // dose distribution image
    int dsize = size[0]*size[1];
    short * dimage = new short[dsize];
    for(int z = 0; z < size[2]; z++) {
      getShortDoseDist(dimage, z);
      ofile.write((char *)dimage, dsize*sizeof(short));

      if(DEBUG || kVerbose > 0) {
	for(int j = 0; j < dsize; j++) {
	  if(dimage[j] < 0)
	    G4cout << "[" << j << "," << z << "]"
		      << dimage[j] << ", ";
	}
      }
    }
    if(DEBUG || kVerbose > 0) G4cout << G4endl;
    delete [] dimage;

    // relative location of the dose distribution image for 
    // the modality image
    kDose[0].getCenterPosition(fCenter);
    for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
    ofile.write((char *)iCenter, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image relative location : ("
		<< iCenter[0] << ", "
		<< iCenter[1] << ", "
		<< iCenter[2] << ")" << G4endl;
    }

  }

  //----- ROI image -----//
  if(!isROIEmpty()) {
    // ROI image size
    kRoi[0].getSize(size);
    ofile.write((char *)size, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    kRoi[0].getMinMax(minmax);
    ofile.write((char *)minmax, sizeof(short)*2);

    // ROI distribution scaling 
    scale = (float)kRoi[0].getScale();
    ofile.write((char *)&scale, sizeof(float));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    for(int i = 0; i < size[2]; i++) {
      short * rimage = kRoi[0].getImage(i);
      ofile.write((char *)rimage, rsize*sizeof(short));

    }

    // ROI relative location
    kRoi[0].getCenterPosition(fCenter);
    for(int i = 0; i < 3; i++) iCenter[i] = (int)fCenter[i];
    ofile.write((char *)iCenter, 3*sizeof(int));
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< iCenter[0] << ", "
		<< iCenter[1] << ", "
		<< iCenter[2] << ")" << G4endl;
    }
  }


  //----- track information -----//
  // track
  int ntrk = (int)kSteps.size();
  ofile.write((char *)&ntrk, sizeof(int));
  if(DEBUG || kVerbose > 0) {
    G4cout << "# of tracks : "
	      << ntrk << G4endl;
  }
  for(int i = 0; i < ntrk; i++) {
    float * tp = kSteps[i];
    ofile.write((char *)tp, sizeof(float)*6);
  }


  // file end mark
  ofile.write("END", 3);

  ofile.close();

  return true;
}
//
bool G4GMocrenIO::storeData2(char * _filename) {
  kFileName = _filename;
  return storeData();
}

bool G4GMocrenIO::retrieveData() {

  // input file open
  std::ifstream ifile(kFileName.c_str(), std::ios_base::in|std::ios_base::binary);
  if(!ifile) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Cannot open file: " << kFileName
	     << " in G4GMocrenIO::retrieveData()." << G4endl;
    return false;
  }

  // file identifier
  char verid[9];
  ifile.read((char *)verid, 8);
  // file version
  unsigned char ver;
  ifile.read((char *)&ver, 1);
  ifile.close();

  if(std::strncmp(verid, "gMocren", 7) == 0) {
    if(ver == 0x03) {
      G4cout << ">>>>>>>  retrieve data (ver.3) <<<<<<<" << G4endl;
      G4cout << "         " << kFileName << G4endl;
      retrieveData3();
    } else if (ver == 0x04) {
      G4cout << ">>>>>>>  retrieve data (ver.4) <<<<<<<" << G4endl;
      G4cout << "         " << kFileName << G4endl;
      retrieveData4();
    } else {
      if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	G4cout << "Error -- invalid file version : " << (int)ver
		  << G4endl;
	G4cout << "         " << kFileName << G4endl;
      }
      G4Exception("G4GMocrenIO::retrieveDadta()",
                  "gMocren2001", FatalException,
                  "Error.");

    }
  } else if(std::strncmp(verid, "GRAPE", 5) == 0) {
    G4cout << ">>>>>>>  retrieve data (ver.2) <<<<<<<" << G4endl;
    G4cout << "         " << kFileName << G4endl;
    retrieveData2();
  } else {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << kFileName << " was not gdd file." << G4endl;
    return false;
  }

  return true;
}

bool G4GMocrenIO::retrieveData(char * _filename) {
  kFileName = _filename;
  return retrieveData();
}

// 
bool G4GMocrenIO::retrieveData4() {

  bool DEBUG = false;//

  // input file open
  std::ifstream ifile(kFileName.c_str(), std::ios_base::in|std::ios_base::binary);
  if(!ifile) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Cannot open file: " << kFileName
		<< " in G4GMocrenIO::retrieveData3()." << G4endl;
    return false;
  }

  // data buffer
  char ctmp[24];

  // file identifier
  char verid[9];
  ifile.read((char *)verid, 8);

  // file version
  unsigned char ver;
  ifile.read((char *)&ver, 1);
  std::stringstream ss;
  ss << (int)ver;
  kVersion = ss.str();
  if(DEBUG || kVerbose > 0) G4cout << "File version : " << kVersion << G4endl;

  // endian
  ifile.read((char *)&kLittleEndianInput, sizeof(char));
  if(DEBUG || kVerbose > 0) {
    G4cout << "Endian : ";
    if(kLittleEndianInput == 1) 
      G4cout << " little" << G4endl;
    else {
      G4cout << " big" << G4endl;
    }
  }

  // comment length (fixed size)
  int clength;
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, clength);
  // comment
  char cmt[1025];
  ifile.read((char *)cmt, clength);
  std::string scmt = cmt;
  scmt += '\0';
  setComment(scmt);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Data comment : "
	      << kComment << G4endl;
  }

  // voxel spacings for all images
  ifile.read((char *)ctmp, 12);
  convertEndian(ctmp, kVoxelSpacing[0]);
  convertEndian(ctmp+4, kVoxelSpacing[1]);
  convertEndian(ctmp+8, kVoxelSpacing[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }


  // offset from file starting point to the modality image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToModalityData);

  // # of dose distributions
  ifile.read((char *)ctmp, 4);
  int nDoseDist;
  convertEndian(ctmp, nDoseDist);
  
  // offset from file starting point to the dose image data
  for(int i = 0; i < nDoseDist; i++) {
    ifile.read((char *)ctmp, 4);
    unsigned int dptr;
    convertEndian(ctmp, dptr);
    addPointerToDoseDistData(dptr);
  }

  // offset from file starting point to the ROI image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToROIData);

  // offset from file starting point to the track data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToTrackData);

  // offset from file starting point to the detector data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToDetectorData);

  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", ";
    for(int i = 0; i < nDoseDist; i++)
      G4cout << kPointerToDoseDistData[i] << ", ";
    G4cout << kPointerToROIData << ", "
	      << kPointerToTrackData << ", "
	      << kPointerToDetectorData
	      << G4endl;
  }



  if(kPointerToModalityData == 0 && kPointerToDoseDistData.size() == 0 &&
     kPointerToROIData == 0 && kPointerToTrackData == 0) {
    if(DEBUG || kVerbose > 0) {
      G4cout << "No data." << G4endl;
    }
    return false;
  }

  // event number
  /* ver 1
     ifile.read(ctmp, sizeof(int));
     convertEndian(ctmp, numberOfEvents);
  */

  int size[3];
  float scale;
  double dscale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];

  //----- Modality image -----//
  // modality image size
  ifile.read(ctmp, 3*sizeof(int));
  convertEndian(ctmp, size[0]);
  convertEndian(ctmp+sizeof(int), size[1]);
  convertEndian(ctmp+2*sizeof(int), size[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }
  kModality.setSize(size);

  // modality image voxel spacing
  /*
    ifile.read(ctmp, 3*sizeof(float));
    convertEndian(ctmp, modalityImageVoxelSpacing[0]);
    convertEndian(ctmp+sizeof(float), modalityImageVoxelSpacing[1]);
    convertEndian(ctmp+2*sizeof(float), modalityImageVoxelSpacing[2]);
  */

  if(kPointerToModalityData != 0) {

    // modality density max. & min.
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);
    kModality.setMinMax(minmax);

    // modality image unit
    char munit[13];
    munit[12] = '\0';
    ifile.read((char *)munit, 12);
    std::string smunit = munit;
    setModalityImageUnit(smunit);

    // modality density scale
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, scale);
    kModality.setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Modality image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // modality density
    int psize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Modality image (" << psize << "): ";
    char * cimage = new char[psize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)cimage, psize*sizeof(short));
      short * mimage = new short[psize];
      for(int j = 0; j < psize; j++) {
	convertEndian(cimage+j*sizeof(short), mimage[j]);
      }
      kModality.addImage(mimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << mimage[(size_t)(psize*0.55)] << ", ";
    }
    if(DEBUG || kVerbose > 0) G4cout << G4endl;
    delete [] cimage;

    // modality desity map for CT value
    size_t msize = minmax[1]-minmax[0]+1;
    if(DEBUG || kVerbose > 0) G4cout << "msize: " << msize << G4endl;
    char * pdmap = new char[msize*sizeof(float)];
    ifile.read((char *)pdmap, msize*sizeof(float));
    float ftmp;
    for(int i = 0; i < (int)msize; i++) {
      convertEndian(pdmap+i*sizeof(float), ftmp);
      kModalityImageDensityMap.push_back(ftmp); 
    }
    delete [] pdmap;

    if(DEBUG || kVerbose > 0) {
      G4cout << "density map : " << std::ends;
      for(int i = 0; i < 10; i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
      for(int i = 0; i < 10; i++) G4cout << "..";
      G4cout << G4endl;
      for(size_t i =kModalityImageDensityMap.size() - 10; i <kModalityImageDensityMap.size(); i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
    }

  }


  //----- dose distribution image -----//
  for(int ndose = 0; ndose < nDoseDist; ndose++) {

    newDoseDist();

    // dose distrbution image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }
    kDose[ndose].setSize(size);

    // dose distribution max. & min. 
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);

    // dose distribution unit
    char dunit[13];
    dunit[12] = '\0';
    ifile.read((char *)dunit, 12);
    std::string sdunit = dunit;
    setDoseDistUnit(sdunit, ndose);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. unit : " << kDoseUnit << G4endl;
    }

    // dose distribution scaling 
    ifile.read((char *)ctmp, 4); // sizeof(float)
    convertEndian(ctmp, scale);
    kDose[ndose].setScale(dscale = scale);

    double dminmax[2];
    for(int i = 0; i < 2; i++) dminmax[i] = minmax[i]*dscale;
    kDose[ndose].setMinMax(dminmax);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image min., max., scale : "
		<< dminmax[0] << ", "
		<< dminmax[1] << ", "
		<< scale << G4endl;
    }

    // dose distribution image
    int dsize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Dose dist. (" << dsize << "): ";
    char * di = new char[dsize*sizeof(short)];
    short * shimage = new short[dsize];
    for(int z = 0; z < size[2]; z++) {
      ifile.read((char *)di, dsize*sizeof(short));
      double * dimage = new double[dsize];
      for(int xy = 0; xy < dsize; xy++) {
	convertEndian(di+xy*sizeof(short), shimage[xy]);
	dimage[xy] = shimage[xy]*dscale;
      }
      kDose[ndose].addImage(dimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << z << "]" << dimage[(size_t)(dsize*0.55)] << ", ";

      if(DEBUG || kVerbose > 0) {
	for(int j = 0; j < dsize; j++) {
	  if(dimage[j] < 0)
	    G4cout << "[" << j << "," << z << "]"
		      << dimage[j] << ", ";
	}
      }
    }
    delete [] shimage;
    delete [] di;
    if(DEBUG || kVerbose > 0) G4cout << G4endl;

    ifile.read((char *)ctmp, 3*4); // 3*sizeof(int)
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+4, iCenter[1]);
    convertEndian(ctmp+8, iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = (float)iCenter[i];
    kDose[ndose].setCenterPosition(fCenter);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }


    // dose distribution name
    char cname[81];
    ifile.read((char *)cname, 80);
    std::string dosename = cname;
    setDoseDistName(dosename, ndose);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. name : " << dosename << G4endl;
    }

  }

  //----- ROI image -----//
  if(kPointerToROIData != 0) {

    newROI();

    // ROI image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    kRoi[0].setSize(size);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+sizeof(short), minmax[1]);
    kRoi[0].setMinMax(minmax);

    // ROI distribution scaling 
    ifile.read((char *)ctmp, sizeof(float));
    convertEndian(ctmp, scale);
    kRoi[0].setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    char * ri = new char[rsize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)ri, rsize*sizeof(short));
      short * rimage = new short[rsize];
      for(int j = 0; j < rsize; j++) {
	convertEndian(ri+j*sizeof(short), rimage[j]);
      }
      kRoi[0].addImage(rimage);

    }
    delete [] ri;

    // ROI relative location
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+sizeof(int), iCenter[1]);
    convertEndian(ctmp+2*sizeof(int), iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = iCenter[i];
    kRoi[0].setCenterPosition(fCenter);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }

  }

  //----- track information -----//
  if(kPointerToTrackData != 0) {

    // track
    ifile.read((char *)ctmp, sizeof(int));
    int ntrk;
    convertEndian(ctmp, ntrk);
    if(DEBUG || kVerbose > 0) {
      G4cout << "# of tracks: " << ntrk << G4endl;
    }

    // track position
    unsigned char rgb[3];
    for(int i = 0; i < ntrk; i++) {


      // # of steps in a track
      ifile.read((char *)ctmp, sizeof(int));
      int nsteps;
      convertEndian(ctmp, nsteps);
      
      // track color
      ifile.read((char *)rgb, 3);

      std::vector<float *> steps;
      // steps
      for(int j = 0; j < nsteps; j++) {

	float * steppoint = new float[6];
	ifile.read((char *)ctmp, sizeof(float)*6);

	for(int k = 0; k < 6; k++) {
	  convertEndian(ctmp+k*sizeof(float), steppoint[k]);
	}
	
	steps.push_back(steppoint);
      }

      // add a track to the track container
      addTrack(steps, rgb);

      if(DEBUG || kVerbose > 0) {
	if(i < 5) {
	  G4cout << i << ": " ;
	  for(int j = 0; j < 3; j++) G4cout << steps[0][j] << " ";
	  int nstp = (int)steps.size();
	  G4cout << "<-> ";
	  for(int j = 3; j < 6; j++) G4cout << steps[nstp-1][j] << " ";
	  G4cout << "    rgb( ";
	  for(int j = 0; j < 3; j++) G4cout << (int)rgb[j] << " ";
	  G4cout << ")" << G4endl;
	}
      }
    }


  }


  //----- detector information -----//
  if(kPointerToDetectorData != 0) {

    // number of detectors
    ifile.read((char *)ctmp, sizeof(int));
    int ndet;
    convertEndian(ctmp, ndet);

    if(DEBUG || kVerbose > 0) {
      G4cout << "# of detectors : "
		<< ndet << G4endl;
    }

    for(int nd = 0; nd < ndet; nd++) {

      // # of edges of a detector
      ifile.read((char *)ctmp, sizeof(int));
      int nedges;
      convertEndian(ctmp, nedges);
      if(DEBUG || kVerbose > 0) {
	G4cout << "# of edges in a detector : " << nedges << G4endl;
      }

      // edges
      std::vector<float *> detector;
      char cftmp[24];
      for(int ne = 0; ne < nedges; ne++) {
      
	ifile.read((char *)cftmp, sizeof(float)*6);
	float * edgePoints = new float[6];
	for(int j = 0; j < 6; j++) convertEndian(&cftmp[sizeof(float)*j], edgePoints[j]);
	detector.push_back(edgePoints);

      }

      if(DEBUG || kVerbose > 0) {
	G4cout << " first edge : (" << detector[0][0] << ", "
		  << detector[0][1] << ", "
		  << detector[0][2] << ") - ("
		  << detector[0][3] << ", "
		  << detector[0][4] << ", "
		  << detector[0][5] << ")" << G4endl;
      }

      // detector color
      unsigned char dcolor[3];
      ifile.read((char *)dcolor, 3);
      if(DEBUG || kVerbose > 0) {
	G4cout << " detector color : rgb("
		  << (int)dcolor[0] << ", "
		  << (int)dcolor[1] << ", "
		  << (int)dcolor[2] << G4endl;
      }


      // detector name
      char cname[80];
      ifile.read((char *)cname, 80);
      std::string dname = cname;
      if(DEBUG || kVerbose > 0) {
	G4cout << " detector name : " << dname << G4endl;
      }


      addDetector(dname, detector, dcolor);

    }
  }


  ifile.close();

  return true;
}
bool G4GMocrenIO::retrieveData4(char * _filename) {
  kFileName = _filename;
  return retrieveData();
}

// 
bool G4GMocrenIO::retrieveData3() {

  bool DEBUG = false;//

  // input file open
  std::ifstream ifile(kFileName.c_str(), std::ios_base::in|std::ios_base::binary);
  if(!ifile) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Cannot open file: " << kFileName
		<< " in G4GMocrenIO::retrieveData3()." << G4endl;
    return false;
  }

  // data buffer
  char ctmp[12];

  // file identifier
  char verid[9];
  ifile.read((char *)verid, 8);

  // file version
  unsigned char ver;
  ifile.read((char *)&ver, 1);
  std::stringstream ss;
  ss << (int)ver;
  kVersion = ss.str();
  if(DEBUG || kVerbose > 0) G4cout << "File version : " << kVersion << G4endl;

  // endian
  ifile.read((char *)&kLittleEndianInput, sizeof(char));
  if(DEBUG || kVerbose > 0) {
    G4cout << "Endian : ";
    if(kLittleEndianInput == 1) 
      G4cout << " little" << G4endl;
    else {
      G4cout << " big" << G4endl;
    }
  }

  // comment length (fixed size)
  int clength;
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, clength);
  // comment
  char cmt[1025];
  ifile.read((char *)cmt, clength);
  std::string scmt = cmt;
  setComment(scmt);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Data comment : "
	      << kComment << G4endl;
  }

  // voxel spacings for all images
  ifile.read((char *)ctmp, 12);
  convertEndian(ctmp, kVoxelSpacing[0]);
  convertEndian(ctmp+4, kVoxelSpacing[1]);
  convertEndian(ctmp+8, kVoxelSpacing[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }


  // offset from file starting point to the modality image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToModalityData);

  // # of dose distributions
  ifile.read((char *)ctmp, 4);
  int nDoseDist;
  convertEndian(ctmp, nDoseDist);
  
  // offset from file starting point to the dose image data
  for(int i = 0; i < nDoseDist; i++) {
    ifile.read((char *)ctmp, 4);
    unsigned int dptr;
    convertEndian(ctmp, dptr);
    addPointerToDoseDistData(dptr);
  }

  // offset from file starting point to the ROI image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToROIData);

  // offset from file starting point to the track data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToTrackData);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", ";
    for(int i = 0; i < nDoseDist; i++)
      G4cout << kPointerToDoseDistData[0] << ", ";
    G4cout << kPointerToROIData << ", "
	      << kPointerToTrackData << G4endl;
  }

  if(kPointerToModalityData == 0 && kPointerToDoseDistData.size() == 0 &&
     kPointerToROIData == 0 && kPointerToTrackData == 0) {
    if(DEBUG || kVerbose > 0) {
      G4cout << "No data." << G4endl;
    }
    return false;
  }

  // event number
  /* ver 1
     ifile.read(ctmp, sizeof(int));
     convertEndian(ctmp, numberOfEvents);
  */

  int size[3];
  float scale;
  double dscale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];

  //----- Modality image -----//
  // modality image size
  ifile.read(ctmp, 3*sizeof(int));
  convertEndian(ctmp, size[0]);
  convertEndian(ctmp+sizeof(int), size[1]);
  convertEndian(ctmp+2*sizeof(int), size[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }
  kModality.setSize(size);

  // modality image voxel spacing
  /*
    ifile.read(ctmp, 3*sizeof(float));
    convertEndian(ctmp, modalityImageVoxelSpacing[0]);
    convertEndian(ctmp+sizeof(float), modalityImageVoxelSpacing[1]);
    convertEndian(ctmp+2*sizeof(float), modalityImageVoxelSpacing[2]);
  */

  if(kPointerToModalityData != 0) {

    // modality density max. & min.
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);
    kModality.setMinMax(minmax);

    // modality image unit
    char munit[13];
    ifile.read((char *)munit, 12);
    std::string smunit = munit;
    setModalityImageUnit(smunit);

    // modality density scale
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, scale);
    kModality.setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Modality image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // modality density
    int psize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Modality image (" << psize << "): ";
    char * cimage = new char[psize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)cimage, psize*sizeof(short));
      short * mimage = new short[psize];
      for(int j = 0; j < psize; j++) {
	convertEndian(cimage+j*sizeof(short), mimage[j]);
      }
      kModality.addImage(mimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << mimage[(size_t)(psize*0.55)] << ", ";
    }
    if(DEBUG || kVerbose > 0) G4cout << G4endl;
    delete [] cimage;

    // modality desity map for CT value
    size_t msize = minmax[1]-minmax[0]+1;
    if(DEBUG || kVerbose > 0) G4cout << "msize: " << msize << G4endl;
    char * pdmap = new char[msize*sizeof(float)];
    ifile.read((char *)pdmap, msize*sizeof(float));
    float ftmp;
    for(int i = 0; i < (int)msize; i++) {
      convertEndian(pdmap+i*sizeof(float), ftmp);
      kModalityImageDensityMap.push_back(ftmp); 
    }
    delete [] pdmap;
    if(DEBUG || kVerbose > 0) {
      G4cout << "density map : " << std::ends;
      for(int i = 0; i < 10; i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
      for(int i = 0; i < 10; i++) G4cout << "..";
      G4cout << G4endl;
      for(size_t i =kModalityImageDensityMap.size() - 10; i <kModalityImageDensityMap.size(); i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
    }

  }


  //----- dose distribution image -----//
  for(int ndose = 0; ndose < nDoseDist; ndose++) {

    newDoseDist();

    // dose distrbution image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }
    kDose[ndose].setSize(size);

    // dose distribution max. & min. 
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);

    // dose distribution unit
    char dunit[13];
    ifile.read((char *)dunit, 12);
    std::string sdunit = dunit;
    setDoseDistUnit(sdunit, ndose);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. unit : " << kDoseUnit << G4endl;
    }

    // dose distribution scaling 
    ifile.read((char *)ctmp, 4); // sizeof(float)
    convertEndian(ctmp, scale);
    kDose[ndose].setScale(dscale = scale);

    double dminmax[2];
    for(int i = 0; i < 2; i++) dminmax[i] = minmax[i]*dscale;
    kDose[ndose].setMinMax(dminmax);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image min., max., scale : "
		<< dminmax[0] << ", "
		<< dminmax[1] << ", "
		<< scale << G4endl;
    }

    // dose distribution image
    int dsize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Dose dist. (" << dsize << "): ";
    char * di = new char[dsize*sizeof(short)];
    short * shimage = new short[dsize];
    for(int z = 0; z < size[2]; z++) {
      ifile.read((char *)di, dsize*sizeof(short));
      double * dimage = new double[dsize];
      for(int xy = 0; xy < dsize; xy++) {
	convertEndian(di+xy*sizeof(short), shimage[xy]);
	dimage[xy] = shimage[xy]*dscale;
      }
      kDose[ndose].addImage(dimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << z << "]" << dimage[(size_t)(dsize*0.55)] << ", ";

      if(DEBUG || kVerbose > 0) {
	for(int j = 0; j < dsize; j++) {
	  if(dimage[j] < 0)
	    G4cout << "[" << j << "," << z << "]"
		      << dimage[j] << ", ";
	}
      }
    }
    delete [] shimage;
    delete [] di;
    if(DEBUG || kVerbose > 0) G4cout << G4endl;

    ifile.read((char *)ctmp, 3*4); // 3*sizeof(int)
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+4, iCenter[1]);
    convertEndian(ctmp+8, iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = (float)iCenter[i];
    kDose[ndose].setCenterPosition(fCenter);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }


  }

  //----- ROI image -----//
  if(kPointerToROIData != 0) {

    newROI();

    // ROI image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    kRoi[0].setSize(size);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+sizeof(short), minmax[1]);
    kRoi[0].setMinMax(minmax);

    // ROI distribution scaling 
    ifile.read((char *)ctmp, sizeof(float));
    convertEndian(ctmp, scale);
    kRoi[0].setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    char * ri = new char[rsize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)ri, rsize*sizeof(short));
      short * rimage = new short[rsize];
      for(int j = 0; j < rsize; j++) {
	convertEndian(ri+j*sizeof(short), rimage[j]);
      }
      kRoi[0].addImage(rimage);

    }
    delete [] ri;

    // ROI relative location
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+sizeof(int), iCenter[1]);
    convertEndian(ctmp+2*sizeof(int), iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = iCenter[i];
    kRoi[0].setCenterPosition(fCenter);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }

  }

  //----- track information -----//
  if(kPointerToTrackData != 0) {

    // track
    ifile.read((char *)ctmp, sizeof(int));
    int ntrk;
    convertEndian(ctmp, ntrk);
    if(DEBUG || kVerbose > 0) {
      G4cout << "# of tracks: " << ntrk << G4endl;
    }

    // v4
    std::vector<float *> trkv4;

    // track position
    for(int i = 0; i < ntrk; i++) {
      float * tp = new float[6];

      ifile.read((char *)ctmp, sizeof(float)*3);
      if(DEBUG || kVerbose > 0) if(i < 10) G4cout << i << ": " ;
      for(int j = 0; j < 3; j++) {
	convertEndian(ctmp+j*sizeof(float), tp[j]);
	if(DEBUG || kVerbose > 0) if(i < 10) G4cout << tp[j] << ", ";
      }

      ifile.read((char *)ctmp, sizeof(float)*3);
      for(int j = 0; j < 3; j++) {
	convertEndian(ctmp+j*sizeof(float), tp[j+3]);
	if(DEBUG || kVerbose > 0) if(i < 10) G4cout << tp[j+3] << ", ";
      }
      addTrack(tp);
      if(DEBUG || kVerbose > 0) if(i < 10) G4cout << G4endl;

      // v4
      trkv4.push_back(tp);
    }

    //v4
    unsigned char trkcolorv4[3];

    // track color
    for(int i = 0; i < ntrk; i++) {
      unsigned char * rgb = new unsigned char[3];
      ifile.read((char *)rgb, 3);
      addTrackColor(rgb);

      // v4
      for(int j = 0; j < 3; j++) trkcolorv4[j] = rgb[j];
      std::vector<float *> trk;
      trk.push_back(trkv4[i]);
      addTrack(trk, trkcolorv4);

    }

  }

  ifile.close();

  return true;
}
bool G4GMocrenIO::retrieveData3(char * _filename) {
  kFileName = _filename;
  return retrieveData();
}

// 
bool G4GMocrenIO::retrieveData2() {

  bool DEBUG = false;//

  // input file open
  std::ifstream ifile(kFileName.c_str(), std::ios_base::in|std::ios_base::binary);
  if(!ifile) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "Cannot open file: " << kFileName
		<< " in G4GMocrenIO::retrieveData2()." << G4endl;
    return false;
  }

  // data buffer
  char ctmp[12];

  // file identifier
  char verid[9];
  ifile.read((char *)verid, 8);

  // file version
  unsigned char ver;
  ifile.read((char *)&ver, 1);
  std::stringstream ss;
  ss << (int)ver;
  kVersion = ss.str();
  if(DEBUG || kVerbose > 0) G4cout << "File version : " << kVersion << G4endl;

  // id of version 1
  char idtmp[IDLENGTH];
  ifile.read((char *)idtmp, IDLENGTH);
  kId = idtmp;
  // version of version 1
  char vertmp[VERLENGTH];
  ifile.read((char *)vertmp, VERLENGTH);

  // endian
  ifile.read((char *)&kLittleEndianInput, sizeof(char));
  if(DEBUG || kVerbose > 0) {
    G4cout << "Endian : ";
    if(kLittleEndianInput == 1) 
      G4cout << " little" << G4endl;
    else {
      G4cout << " big" << G4endl;
    }
  }

  // voxel spacings for all images
  ifile.read((char *)ctmp, 12);
  convertEndian(ctmp, kVoxelSpacing[0]);
  convertEndian(ctmp+4, kVoxelSpacing[1]);
  convertEndian(ctmp+8, kVoxelSpacing[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Voxel spacing : ("
	      << kVoxelSpacing[0] << ", "
	      << kVoxelSpacing[1] << ", "
	      << kVoxelSpacing[2]
	      << ") mm " << G4endl;
  }


  // offset from file starting point to the modality image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToModalityData);

  // offset from file starting point to the dose image data
  unsigned int ptddd;
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, ptddd);
  kPointerToDoseDistData.push_back(ptddd);

  // offset from file starting point to the ROI image data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToROIData);

  // offset from file starting point to the track data
  ifile.read((char *)ctmp, 4);
  convertEndian(ctmp, kPointerToTrackData);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Each pointer to data : "
	      << kPointerToModalityData << ", "
	      << kPointerToDoseDistData[0] << ", "
	      << kPointerToROIData << ", "
	      << kPointerToTrackData << G4endl;
  }

  if(kPointerToModalityData == 0 && kPointerToDoseDistData.size() == 0 &&
     kPointerToROIData == 0 && kPointerToTrackData == 0) {
    if(DEBUG || kVerbose > 0) {
      G4cout << "No data." << G4endl;
    }
    return false;
  }

  // event number
  /* ver 1
     ifile.read(ctmp, sizeof(int));
     convertEndian(ctmp, numberOfEvents);
  */

  int size[3];
  float scale;
  double dscale;
  short minmax[2];
  float fCenter[3];
  int iCenter[3];

  //----- Modality image -----//
  // modality image size
  ifile.read(ctmp, 3*sizeof(int));
  convertEndian(ctmp, size[0]);
  convertEndian(ctmp+sizeof(int), size[1]);
  convertEndian(ctmp+2*sizeof(int), size[2]);
  if(DEBUG || kVerbose > 0) {
    G4cout << "Modality image size : ("
	      << size[0] << ", "
	      << size[1] << ", "
	      << size[2] << ")"
	      << G4endl;
  }
  kModality.setSize(size);

  // modality image voxel spacing
  /*
    ifile.read(ctmp, 3*sizeof(float));
    convertEndian(ctmp, modalityImageVoxelSpacing[0]);
    convertEndian(ctmp+sizeof(float), modalityImageVoxelSpacing[1]);
    convertEndian(ctmp+2*sizeof(float), modalityImageVoxelSpacing[2]);
  */

  if(kPointerToModalityData != 0) {

    // modality density max. & min.
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);
    kModality.setMinMax(minmax);

    // modality density scale
    ifile.read((char *)ctmp, 4);
    convertEndian(ctmp, scale);
    kModality.setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Modality image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // modality density
    int psize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Modality image (" << psize << "): ";
    char * cimage = new char[psize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)cimage, psize*sizeof(short));
      short * mimage = new short[psize];
      for(int j = 0; j < psize; j++) {
	convertEndian(cimage+j*sizeof(short), mimage[j]);
      }
      kModality.addImage(mimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << i << "]" << mimage[(size_t)(psize*0.55)] << ", ";
    }
    if(DEBUG || kVerbose > 0) G4cout << G4endl;
    delete [] cimage;

    // modality desity map for CT value
    size_t msize = minmax[1]-minmax[0]+1;
    if(DEBUG || kVerbose > 0) G4cout << "msize: " << msize << G4endl;
    char * pdmap = new char[msize*sizeof(float)];
    ifile.read((char *)pdmap, msize*sizeof(float));
    float ftmp;
    for(int i = 0; i < (int)msize; i++) {
      convertEndian(pdmap+i*sizeof(float), ftmp);
      kModalityImageDensityMap.push_back(ftmp); 
    }
    delete [] pdmap;
    if(DEBUG || kVerbose > 0) {
      G4cout << "density map : " << std::ends;
      for(int i = 0; i < 10; i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
      for(int i = 0; i < 10; i++) G4cout << "..";
      G4cout << G4endl;
      for(size_t i =kModalityImageDensityMap.size() - 10; i <kModalityImageDensityMap.size(); i++)
	G4cout <<kModalityImageDensityMap[i] << ", ";
      G4cout << G4endl;
    }

  }


  //----- dose distribution image -----//
  if(kPointerToDoseDistData[0] != 0) {

    newDoseDist();

    // dose distrbution image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }
    kDose[0].setSize(size);

    // dose distribution max. & min. 
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+2, minmax[1]);
    // dose distribution scaling 
    ifile.read((char *)ctmp, sizeof(float));
    convertEndian(ctmp, scale);
    kDose[0].setScale(dscale = scale);

    double dminmax[2];
    for(int i = 0; i < 2; i++) dminmax[i] = minmax[i]*dscale;
    kDose[0].setMinMax(dminmax);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image min., max., scale : "
		<< dminmax[0] << ", "
		<< dminmax[1] << ", "
		<< scale << G4endl;
    }

    // dose distribution image
    int dsize = size[0]*size[1];
    if(DEBUG || kVerbose > 0) G4cout << "Dose dist. (" << dsize << "): ";
    char * di = new char[dsize*sizeof(short)];
    short * shimage = new short[dsize];
    for(int z = 0; z < size[2]; z++) {
      ifile.read((char *)di, dsize*sizeof(short));
      double * dimage = new double[dsize];
      for(int xy = 0; xy < dsize; xy++) {
	convertEndian(di+xy*sizeof(short), shimage[xy]);
	dimage[xy] = shimage[xy]*dscale;
      }
      kDose[0].addImage(dimage);

      if(DEBUG || kVerbose > 0) G4cout << "[" << z << "]" << dimage[(size_t)(dsize*0.55)] << ", ";

      if(DEBUG || kVerbose > 0) {
	for(int j = 0; j < dsize; j++) {
	  if(dimage[j] < 0)
	    G4cout << "[" << j << "," << z << "]"
		      << dimage[j] << ", ";
	}
      }
    }
    delete [] shimage;
    delete [] di;
    if(DEBUG || kVerbose > 0) G4cout << G4endl;

    /* ver 1
       float doseDist;
       int dosePid;
       double * doseData = new double[numDoseImageVoxels];
       for(int i = 0; i < numDose; i++) {
       ifile.read(ctmp, sizeof(int));
       convertEndian(ctmp, dosePid);
       for(int j = 0; j < numDoseImageVoxels; j++) {
       ifile.read(ctmp, sizeof(float));
       convertEndian(ctmp, doseDist);
       doseData[j] = doseDist;
       }
       setDose(dosePid, doseData);
       }
       delete [] doseData;
       if(totalDose == NULL) totalDose = new double[numDoseImageVoxels];
       for(int i = 0; i < numDoseImageVoxels; i++) {
       ifile.read(ctmp, sizeof(float));
       convertEndian(ctmp, doseDist);
       totalDose[i] = doseDist;
       }
    */

    /* ver 1
    // relative location between the two images
    ifile.read(ctmp, 3*sizeof(float));
    convertEndian(ctmp, relativeLocation[0]);
    convertEndian(ctmp+sizeof(float), relativeLocation[1]);
    convertEndian(ctmp+2*sizeof(float), relativeLocation[2]);
    */

    // relative location of the dose distribution image for 
    // the modality image
    //ofile.write((char *)relativeLocation, 3*sizeof(float));
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+sizeof(int), iCenter[1]);
    convertEndian(ctmp+2*sizeof(int), iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = (float)iCenter[i];
    kDose[0].setCenterPosition(fCenter);

    if(DEBUG || kVerbose > 0) {
      G4cout << "Dose dist. image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }


  }

  //----- ROI image -----//
  if(kPointerToROIData != 0) {

    newROI();

    // ROI image size
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, size[0]);
    convertEndian(ctmp+sizeof(int), size[1]);
    convertEndian(ctmp+2*sizeof(int), size[2]);
    kRoi[0].setSize(size);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image size : ("
		<< size[0] << ", "
		<< size[1] << ", "
		<< size[2] << ")"
		<< G4endl;
    }

    // ROI max. & min.
    ifile.read((char *)ctmp, sizeof(short)*2);
    convertEndian(ctmp, minmax[0]);
    convertEndian(ctmp+sizeof(short), minmax[1]);
    kRoi[0].setMinMax(minmax);

    // ROI distribution scaling 
    ifile.read((char *)ctmp, sizeof(float));
    convertEndian(ctmp, scale);
    kRoi[0].setScale(dscale = scale);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image min., max., scale : "
		<< minmax[0] << ", "
		<< minmax[1] << ", "
		<< scale << G4endl;
    }

    // ROI image
    int rsize = size[0]*size[1];
    char * ri = new char[rsize*sizeof(short)];
    for(int i = 0; i < size[2]; i++) {
      ifile.read((char *)ri, rsize*sizeof(short));
      short * rimage = new short[rsize];
      for(int j = 0; j < rsize; j++) {
	convertEndian(ri+j*sizeof(short), rimage[j]);
      }
      kRoi[0].addImage(rimage);

    }
    delete [] ri;

    // ROI relative location
    ifile.read((char *)ctmp, 3*sizeof(int));
    convertEndian(ctmp, iCenter[0]);
    convertEndian(ctmp+sizeof(int), iCenter[1]);
    convertEndian(ctmp+2*sizeof(int), iCenter[2]);
    for(int i = 0; i < 3; i++) fCenter[i] = iCenter[i];
    kRoi[0].setCenterPosition(fCenter);
    if(DEBUG || kVerbose > 0) {
      G4cout << "ROI image relative location : ("
		<< fCenter[0] << ", "
		<< fCenter[1] << ", "
		<< fCenter[2] << ")" << G4endl;
    }

  }

  //----- track information -----//
  if(kPointerToTrackData != 0) {

    // track
    ifile.read((char *)ctmp, sizeof(int));
    int ntrk;
    convertEndian(ctmp, ntrk);
    if(DEBUG || kVerbose > 0) {
      G4cout << "# of tracks: " << ntrk << G4endl;
    }

    //v4
    unsigned char trkcolorv4[3] = {255, 0, 0};

    for(int i = 0; i < ntrk; i++) {
      float * tp = new float[6];
      // v4
      std::vector<float *> trkv4;

      ifile.read((char *)ctmp, sizeof(float)*3);
      if(DEBUG || kVerbose > 0) if(i < 10) G4cout << i << ": " ;
      for(int j = 0; j < 3; j++) {
	convertEndian(ctmp+j*sizeof(float), tp[j]);
	if(DEBUG || kVerbose > 0) if(i < 10) G4cout << tp[j] << ", ";
      }

      ifile.read((char *)ctmp, sizeof(float)*3);
      for(int j = 0; j < 3; j++) {
	convertEndian(ctmp+j*sizeof(float), tp[j+3]);
	if(DEBUG || kVerbose > 0) if(i < 10) G4cout << tp[j+3] << ", ";
      }

      kSteps.push_back(tp);
      // v4
      trkv4.push_back(tp);
      addTrack(trkv4, trkcolorv4);
      
      if(DEBUG || kVerbose > 0) if(i < 10) G4cout << G4endl;
    }

  }

  /* ver 1
  // track
  int ntracks;
  ifile.read(ctmp, sizeof(int));
  convertEndian(ctmp, ntracks);
  // track displacement
  ifile.read(ctmp, 3*sizeof(float));
  convertEndian(ctmp, trackDisplacement[0]);
  convertEndian(ctmp+sizeof(float), trackDisplacement[2]); // exchanged with [1]
  convertEndian(ctmp+2*sizeof(float), trackDisplacement[1]);
  //
  //for(int i = 0; i < ntracks && i < 100; i++) {
  for(int i = 0; i < ntracks; i++) {
  DicomDoseTrack trk;
  short trackid, parentid, pid;
  int npoints;
  ifile.read(ctmp, sizeof(short));
  convertEndian(ctmp, trackid);
  trk.setID(trackid);
  ifile.read(ctmp, sizeof(short));
  convertEndian(ctmp, parentid);
  trk.setParentID(parentid);
  ifile.read(ctmp, sizeof(short));
  convertEndian(ctmp, pid);
  trk.setPID(pid);
  ifile.read(ctmp, sizeof(int));
  convertEndian(ctmp, npoints);
  for(int i = 0; i < npoints; i++) {
  ifile.read(ctmp, 3*sizeof(float));
  // storing only start and end points
  //if(i == 0 || i == npoints - 1) {
  float * point = new float[3];
  convertEndian(ctmp, point[0]);
  convertEndian(ctmp+sizeof(float), point[1]);
  convertEndian(ctmp+2*sizeof(float), point[2]);
  trk.addPoint(point);
  //}
  }
  track.push_back(trk);
  }
  */

  ifile.close();

  return true;
}

bool G4GMocrenIO::retrieveData2(char * _filename) {
  kFileName = _filename;
  return retrieveData();
}

void G4GMocrenIO::setID() {
  time_t t;
  time(&t);

  tm * ti;
  ti = localtime(&t);

  char cmonth[12][4] = {"Jan", "Feb", "Mar", "Apr",
			"May", "Jun", "Jul", "Aug",
			"Sep", "Oct", "Nov", "Dec"};
  std::stringstream ss;
  ss << std::setfill('0')
     << std::setw(2)
     << ti->tm_hour << ":"
     << std::setw(2)
     << ti->tm_min << ":"
     << std::setw(2)
     << ti->tm_sec << ","
     << cmonth[ti->tm_mon] << "."
     << std::setw(2)
     << ti->tm_mday << ","
     << ti->tm_year+1900;

  kId = ss.str();
}

// get & set the file version
std::string & G4GMocrenIO::getVersion() {return kVersion;}
void G4GMocrenIO::setVersion(std::string & _version) {kVersion = _version;}

// set endians of input/output data
void G4GMocrenIO::setLittleEndianInput(bool _little) {kLittleEndianInput = _little;}
void G4GMocrenIO::setLittleEndianOutput(bool _little) {kLittleEndianOutput = _little;}

// voxel spacing
void G4GMocrenIO::setVoxelSpacing(float _spacing[3]) {
  for(int i = 0; i < 3; i++) kVoxelSpacing[i] = _spacing[i];
}
void G4GMocrenIO::getVoxelSpacing(float _spacing[3]) {
  for(int i = 0; i < 3; i++) _spacing[i] = kVoxelSpacing[i];
}

// get & set number of events
int & G4GMocrenIO::getNumberOfEvents() {
  return kNumberOfEvents;
}
void G4GMocrenIO::setNumberOfEvents(int & _numberOfEvents) {
  kNumberOfEvents = _numberOfEvents;
}
void G4GMocrenIO::addOneEvent() {
  kNumberOfEvents++;
}

// set/get pointer the modality image data
void G4GMocrenIO::setPointerToModalityData(unsigned int & _pointer) {
  kPointerToModalityData = _pointer;
}
unsigned int G4GMocrenIO::getPointerToModalityData() {
  return kPointerToModalityData;
}
// set/get pointer the dose distribution image data
void G4GMocrenIO::addPointerToDoseDistData(unsigned int & _pointer) {
  kPointerToDoseDistData.push_back(_pointer);
}
unsigned int G4GMocrenIO::getPointerToDoseDistData(int _elem) {
  if(kPointerToDoseDistData.size() == 0 ||
     kPointerToDoseDistData.size() < (size_t)_elem)
    return 0;
  else
    return kPointerToDoseDistData[_elem];
}

// set/get pointer the ROI image data
void G4GMocrenIO::setPointerToROIData(unsigned int & _pointer) {
  kPointerToROIData = _pointer;
}
unsigned int G4GMocrenIO::getPointerToROIData() {
  return kPointerToROIData;
}
// set/get pointer the track data
void G4GMocrenIO::setPointerToTrackData(unsigned int & _pointer) {
  kPointerToTrackData = _pointer;
}
unsigned int G4GMocrenIO::getPointerToTrackData() {
  return kPointerToTrackData;
}

// calculate pointers for version 4
void G4GMocrenIO::calcPointers4() {

  // pointer to modality data
  unsigned int pointer = 1070; // up to "pointer to the detector data" except for "pointer to the dose dist data"
  int nDoseDist = getNumDoseDist();
  pointer += nDoseDist*4;

  setPointerToModalityData(pointer);

  // pointer to dose data
  // ct-density map for modality data
  int msize[3];
  getModalityImageSize(msize);
  short mminmax[2];
  getModalityImageMinMax(mminmax);
  int pmsize = 2*msize[0]*msize[1]*msize[2];
  int pmmap = 4*(mminmax[1] - mminmax[0] + 1);
  pointer += 32 + pmsize + pmmap;
  //
  kPointerToDoseDistData.clear();
  if(nDoseDist == 0) {
    unsigned int pointer0 = 0;
    addPointerToDoseDistData(pointer0);
  }
  for(int ndose = 0; ndose < nDoseDist; ndose++) {
    addPointerToDoseDistData(pointer);
    int dsize[3];
    getDoseDistSize(dsize);
    pointer += 44 + dsize[0]*dsize[1]*dsize[2]*2 + 80;
  }

  // pointer to roi data
  if(!isROIEmpty()) {
    setPointerToROIData(pointer);
    
    int rsize[3];
    getROISize(rsize);
    int prsize = 2*rsize[0]*rsize[1]*rsize[2];
    pointer += 20 + prsize + 12;
  } else {
    unsigned int pointer0 = 0;
    setPointerToROIData(pointer0);
  }

  // pointer to track data
  int ntrk = (int)kTracks.size();
  if(ntrk != 0) {
    setPointerToTrackData(pointer);

    pointer += 4; // # of tracks
    for(int nt = 0; nt < ntrk; nt++) {
      int nsteps = kTracks[nt].getNumberOfSteps();
      pointer += 4 + 3 + nsteps*(4*6); // # of steps + color + steps(float*6)
    }
  } else {
    unsigned int pointer0 = 0;
    setPointerToTrackData(pointer0);
  }
  if(kVerbose > 0) G4cout << " pointer to the track data :"
			     << kPointerToTrackData << G4endl;

  // pointer to detector data
  int ndet = (int)kDetectors.size();
  if(ndet != 0) {
    kPointerToDetectorData = pointer;
  } else {
    kPointerToDetectorData = 0;
  }
  if(kVerbose > 0) G4cout << " pointer to the detector data :"
			     << kPointerToDetectorData << G4endl;

}

// calculate pointers for ver.3
void G4GMocrenIO::calcPointers3() {

  // pointer to modality data
  unsigned int pointer = 1066; // up to "pointer to the track data" except for "pointer to the dose dist data"
  int nDoseDist = getNumDoseDist();
  pointer += nDoseDist*4;

  setPointerToModalityData(pointer);

  // pointer to dose data
  // ct-density map for modality data
  int msize[3];
  getModalityImageSize(msize);
  short mminmax[2];
  getModalityImageMinMax(mminmax);
  int pmsize = 2*msize[0]*msize[1]*msize[2];
  int pmmap = 4*(mminmax[1] - mminmax[0] + 1);
  pointer += 32 + pmsize + pmmap;
  //
  kPointerToDoseDistData.clear();
  if(nDoseDist == 0) {
    unsigned int pointer0 = 0;
    addPointerToDoseDistData(pointer0);
  }
  for(int ndose = 0; ndose < nDoseDist; ndose++) {
    addPointerToDoseDistData(pointer);
    int dsize[3];
    getDoseDistSize(dsize);
    pointer += 44 + dsize[0]*dsize[1]*dsize[2]*2;
  }

  // pointer to roi data
  if(!isROIEmpty()) {
    setPointerToROIData(pointer);
    
    int rsize[3];
    getROISize(rsize);
    int prsize = 2*rsize[0]*rsize[1]*rsize[2];
    pointer += 20 + prsize + 12;
  } else {
    unsigned int pointer0 = 0;
    setPointerToROIData(pointer0);
  }

  //
  if(getNumTracks() != 0) 
    setPointerToTrackData(pointer);
  else {
    unsigned int pointer0 = 0;
    setPointerToTrackData(pointer0);
  }

}

// calculate pointers for ver.2
void G4GMocrenIO::calcPointers2() {

  // pointer to modality data
  unsigned int pointer = 65;
  setPointerToModalityData(pointer);

  // pointer to dose data
  int msize[3];
  getModalityImageSize(msize);
  short mminmax[2];
  getModalityImageMinMax(mminmax);
  int pmsize = 2*msize[0]*msize[1]*msize[2];
  int pmmap = 4*(mminmax[1] - mminmax[0] + 1);
  pointer += 20 + pmsize + pmmap;
  int dsize[3];
  getDoseDistSize(dsize);
  kPointerToDoseDistData.clear();
  if(dsize[0] != 0) {
    kPointerToDoseDistData.push_back(pointer);

    int pdsize = 2*dsize[0]*dsize[1]*dsize[2];
    pointer += 20 + pdsize + 12;
  } else {
    unsigned int pointer0 = 0;
    kPointerToDoseDistData.push_back(pointer0);
  }

  // pointer to roi data
  if(!isROIEmpty())  {
    int rsize[3];
    getROISize(rsize);
    setPointerToROIData(pointer);
    int prsize = 2*rsize[0]*rsize[1]*rsize[2];
    pointer += 20 + prsize + 12;

  } else {
    unsigned int pointer0 = 0;
    setPointerToROIData(pointer0);
  }

  //
  if(getNumTracks() != 0) 
    setPointerToTrackData(pointer);
  else {
    unsigned int pointer0 = 0;
    setPointerToTrackData(pointer0);
  }

}


//----- Modality image -----//
void G4GMocrenIO::getModalityImageSize(int _size[3]) {

  kModality.getSize(_size);
}
void G4GMocrenIO::setModalityImageSize(int _size[3]) {

  kModality.setSize(_size);
}

// get & set the modality image size
void G4GMocrenIO::setModalityImageScale(double & _scale) {

  kModality.setScale(_scale);
}
double G4GMocrenIO::getModalityImageScale() {

  return kModality.getScale();
}

// set the modality image in CT 
void G4GMocrenIO::setModalityImage(short * _image) {

  kModality.addImage(_image);
}
short * G4GMocrenIO::getModalityImage(int _z) {
  
  return kModality.getImage(_z);
}
void G4GMocrenIO::clearModalityImage() {
  
  kModality.clearImage();
}
// set/get the modality image density map
void G4GMocrenIO::setModalityImageDensityMap(std::vector<float> & _map) {
  kModalityImageDensityMap = _map;
}
std::vector<float> & G4GMocrenIO::getModalityImageDensityMap() {
  return kModalityImageDensityMap;
}
// set the modality image min./max.
void G4GMocrenIO::setModalityImageMinMax(short _minmax[2]) {

  kModality.setMinMax(_minmax);
}  
// get the modality image min./max.
void G4GMocrenIO::getModalityImageMinMax(short _minmax[2]) {

  short minmax[2];
  kModality.getMinMax(minmax);
  for(int i = 0; i < 2; i++) _minmax[i] = minmax[i];
}  
short G4GMocrenIO::getModalityImageMax() {

  short minmax[2];
  kModality.getMinMax(minmax);
  return minmax[1];
}
short G4GMocrenIO::getModalityImageMin() {

  short minmax[2];
  kModality.getMinMax(minmax);
  return minmax[0];
}
// set/get position of the modality image center
void G4GMocrenIO::setModalityCenterPosition(float _center[3]) {

  kModality.setCenterPosition(_center);
}
void G4GMocrenIO::getModalityCenterPosition(float _center[3]) {

  if(isROIEmpty())
    for(int i = 0; i < 3; i++) _center[i] = 0;
  else 
    kModality.getCenterPosition(_center);
}
// get & set the modality image unit
std::string G4GMocrenIO::getModalityImageUnit() {
  return kModalityUnit;
}
void G4GMocrenIO::setModalityImageUnit(std::string & _unit) {
  kModalityUnit = _unit;
}
//
short G4GMocrenIO::convertDensityToHU(float & _dens) {
  short rval = -1024; // default: air
  int nmap = (int)kModalityImageDensityMap.size();
  if(nmap != 0) {
    short minmax[2];
    kModality.getMinMax(minmax);
    rval = minmax[1];
    for(int i = 0; i < nmap; i++) {
      //G4cout << kModalityImageDensityMap[i] << G4endl;
      if(_dens <= kModalityImageDensityMap[i]) {
	rval = i + minmax[0];
	break;
      }
    }
  }
  return rval;
}


//----- Dose distribution -----//
//
void G4GMocrenIO::newDoseDist() {
  GMocrenDataPrimitive<double> doseData;
  kDose.push_back(doseData);
}
int G4GMocrenIO::getNumDoseDist() {
  return (int)kDose.size();
}

// get & set the dose distribution unit
std::string G4GMocrenIO::getDoseDistUnit(int _num) {
  // to avoid a warning in the compile process
  if(kDoseUnit.size() > static_cast<size_t>(_num)) return kDoseUnit;

  return kDoseUnit;
}
void G4GMocrenIO::setDoseDistUnit(std::string & _unit, int _num) {
  // to avoid a warning in the compile process
  if(_unit.size() > static_cast<size_t>(_num)) kDoseUnit = _unit;

  //char unit[13];
  //std::strncpy(unit, _unit.c_str(), 12);
  //doseUnit = unit;
  kDoseUnit = _unit;
}
//
void G4GMocrenIO::getDoseDistSize(int _size[3], int _num) {
  if(isDoseEmpty())
    for(int i = 0; i < 3; i++) _size[i] = 0;
  else 
    kDose[_num].getSize(_size);
}
void G4GMocrenIO::setDoseDistSize(int _size[3], int _num) {

  kDose[_num].setSize(_size);

  //resetDose();
}

void G4GMocrenIO::setDoseDistMinMax(short _minmax[2], int _num) {

  double minmax[2];
  double scale = kDose[_num].getScale();
  for(int i = 0; i < 2; i++) minmax[i] = (double)_minmax[i]*scale;
  kDose[_num].setMinMax(minmax);
}  
void G4GMocrenIO::getDoseDistMinMax(short _minmax[2], int _num) {

  if(isDoseEmpty())
    for(int i = 0; i < 2; i++) _minmax[i] = 0;
  else {
    double minmax[2];
    kDose[_num].getMinMax(minmax);
    double scale = kDose[_num].getScale();
    for(int i = 0; i < 2; i++) _minmax[i] = (short)(minmax[i]/scale+0.5);
  }
}  
void G4GMocrenIO::setDoseDistMinMax(double _minmax[2], int _num) {

  kDose[_num].setMinMax(_minmax);
}  
void G4GMocrenIO::getDoseDistMinMax(double _minmax[2], int _num) {

  if(isDoseEmpty())
    for(int i = 0; i < 2; i++) _minmax[i] = 0.;
  else
    kDose[_num].getMinMax(_minmax);
}  

// get & set the dose distribution image scale
void G4GMocrenIO::setDoseDistScale(double & _scale, int _num) {

  kDose[_num].setScale(_scale);
}
double G4GMocrenIO::getDoseDistScale(int _num) {

  if(isDoseEmpty())
    return 0.;
  else 
    return kDose[_num].getScale();
}

/*
  void G4GMocrenIO::initializeShortDoseDist() {
  ;
  }
  void G4GMocrenIO::finalizeShortDoseDist() {
  ;
  }
*/
// set the dose distribution image
void G4GMocrenIO::setShortDoseDist(short * _image, int _num) {

  int size[3];
  kDose[_num].getSize(size);
  int dsize = size[0]*size[1];
  double * ddata = new double[dsize];
  double scale = kDose[_num].getScale();
  double minmax[2];
  kDose[_num].getMinMax(minmax);
  for(int xy = 0; xy < dsize; xy++) {
    ddata[xy] = _image[xy]*scale;
    if(ddata[xy] < minmax[0]) minmax[0] = ddata[xy];
    if(ddata[xy] > minmax[1]) minmax[1] = ddata[xy];
  }
  kDose[_num].addImage(ddata);

  // set min./max.
  kDose[_num].setMinMax(minmax);
}
void G4GMocrenIO::getShortDoseDist(short * _data, int _z, int _num) {

  if(_data == NULL) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "In G4GMocrenIO::getShortDoseDist(), "
		<< "first argument is NULL pointer. "
		<< "The argument must be allocated array."
		<< G4endl;
    G4Exception("G4GMocrenIO::getShortDoseDist()",
                "gMocren2002", FatalException,
                "Error.");
    return;
  }

  int size[3];
  kDose[_num].getSize(size);
  //short * shdata = new short[size[0]*size[1]];
  double * ddata = kDose[_num].getImage(_z);
  double scale = kDose[_num].getScale();
  for(int xy = 0; xy < size[0]*size[1]; xy++) {
    _data[xy] = (short)(ddata[xy]/scale+0.5); //there is never negative value
  }
}
void G4GMocrenIO::getShortDoseDistMinMax(short _minmax[2], int _num) {
  double scale = kDose[_num].getScale();
  double minmax[2];
  kDose[_num].getMinMax(minmax);
  for(int i = 0; i < 2; i++)
    _minmax[i] = (short)(minmax[i]/scale+0.5);
}
//
void G4GMocrenIO::setDoseDist(double * _image, int _num) {

  kDose[_num].addImage(_image);
}
double * G4GMocrenIO::getDoseDist(int _z, int _num) {

  double * image;
  if(isDoseEmpty()) {
    image = 0;
  } else {
    image = kDose[_num].getImage(_z);
  }
  return image;
}
/*
  void G4GMocrenIO::getDoseDist(double * & _image, int _z, int _num) {

  G4cout << " <" << (void*)_image << "> ";
  if(isDoseEmpty()) {
  _image = 0;
  } else {
  _image = kDose[_num].getImage(_z);
  G4cout << " <" << (void*)_image << "> ";
  G4cout << _image[100] << " ";
  }
  }
*/
bool G4GMocrenIO::addDoseDist(std::vector<double *> & _image, int _num) {

  int size[3];
  getDoseDistSize(size, _num);
  std::vector<double *> dosedist = kDose[_num].getImage();

  int nimg = size[0]*size[1];
  for(int z = 0; z < size[2]; z++) {
    for(int xy = 0; xy < nimg; xy++) {
      dosedist[z][xy] += _image[z][xy];
    }
  }

  return true;
}
//void setDoseDistDensityMap(float * _map) {doseImageDensityMap = _map;};
// set the dose distribution image displacement
void G4GMocrenIO::setDoseDistCenterPosition(float _center[3], int _num) {

  kDose[_num].setCenterPosition(_center);
}
void G4GMocrenIO::getDoseDistCenterPosition(float _center[3], int _num) {

  if(isDoseEmpty())
    for(int i = 0; i < 3; i++) _center[i] = 0;
  else 
    kDose[_num].getCenterPosition(_center);
}
// set & get name of dose distribution
void G4GMocrenIO::setDoseDistName(std::string _name, int _num) {

  kDose[_num].setName(_name);
}
std::string G4GMocrenIO::getDoseDistName(int _num) {

  std::string name;
  if(isDoseEmpty())
    return name;
  else 
    return kDose[_num].getName();
}
// copy dose distributions
void G4GMocrenIO::copyDoseDist(std::vector<class GMocrenDataPrimitive<double> > & _dose) {
  std::vector<class GMocrenDataPrimitive<double> >::iterator itr;
  for(itr = kDose.begin(); itr != kDose.end(); itr++) {
    _dose.push_back(*itr);
  }
}
// merge two dose distributions
bool G4GMocrenIO::mergeDoseDist(std::vector<class GMocrenDataPrimitive<double> > & _dose) {
  if(kDose.size() != _dose.size()) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
      G4cout << "G4GMocrenIO::mergeDoseDist() : Error" << G4endl; 
      G4cout << "   Unable to merge the dose distributions,"<< G4endl;
      G4cout << "   because of different size of dose maps."<< G4endl;
    }
    return false;
  }

  int num = (int)kDose.size();
  std::vector<class GMocrenDataPrimitive<double> >::iterator itr1 = kDose.begin();
  std::vector<class GMocrenDataPrimitive<double> >::iterator itr2 = _dose.begin();
  for(int i = 0; i < num; i++, itr1++, itr2++) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      if(kVerbose > 0)
	G4cout << "merged dose distribution [" << i << "]" << G4endl;
    *itr1 += *itr2;
  }

  return true;
}
//
void G4GMocrenIO::clearDoseDistAll() {

  if(!isDoseEmpty()) {
    for(int i = 0; i < getNumDoseDist(); i++) {
      kDose[i].clear();
    }
    kDose.clear();
  }
}
//
bool G4GMocrenIO::isDoseEmpty() {
  if(kDose.empty()) {
    //if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    //  G4cout << "!!! dose distribution data is empty." << G4endl;
    return true;
  } else {
    return false;
  }
}

//
void G4GMocrenIO::calcDoseDistScale() {

  double scale;
  double minmax[2];

  for(int i = 0; i < (int)kDose.size(); i++) {
    kDose[i].getMinMax(minmax);
    scale = minmax[1]/DOSERANGE;
    kDose[i].setScale(scale);
  }
}


//----- RoI -----//

// add one RoI data
void G4GMocrenIO::newROI() {
  GMocrenDataPrimitive<short>  roiData;
  kRoi.push_back(roiData);
}
int G4GMocrenIO::getNumROI() {
  return (int)kRoi.size();
}

// set/get the ROI image scale
void G4GMocrenIO::setROIScale(double & _scale, int _num) {

  kRoi[_num].setScale(_scale);
}
double G4GMocrenIO::getROIScale(int _num) {

  if(isROIEmpty())
    return 0.;
  else 
    return kRoi[_num].getScale();
}
// set the ROI image 
void G4GMocrenIO::setROI(short * _image, int _num) {

  kRoi[_num].addImage(_image);
}
short * G4GMocrenIO::getROI(int _z, int _num) {

  if(isROIEmpty())
    return 0;
  else 
    return kRoi[_num].getImage(_z);
}
// set/get the ROI image size
void G4GMocrenIO::setROISize(int _size[3], int _num) {

  return kRoi[_num].setSize(_size);
}
void G4GMocrenIO::getROISize(int _size[3], int _num) {

  if(isROIEmpty())
    for(int i = 0; i < 3; i++) _size[i] = 0;
  else 
    return kRoi[_num].getSize(_size);
}
// set/get the ROI image min. and max.
void G4GMocrenIO::setROIMinMax(short _minmax[2], int _num) {

  kRoi[_num].setMinMax(_minmax);
}
void G4GMocrenIO::getROIMinMax(short _minmax[2], int _num) {

  if(isROIEmpty())
    for(int i = 0; i < 2; i++) _minmax[i] = 0;
  else 
    kRoi[_num].getMinMax(_minmax);
}
// set/get the ROI image displacement
void G4GMocrenIO::setROICenterPosition(float _center[3], int _num) {

  kRoi[_num].setCenterPosition(_center);
}
void G4GMocrenIO::getROICenterPosition(float _center[3], int _num) {

  if(isROIEmpty())
    for(int i = 0; i < 3; i++) _center[i] = 0;
  else 
    kRoi[_num].getCenterPosition(_center);
}
//
void G4GMocrenIO::clearROIAll() {

  if(!isROIEmpty()) {
    for(int i = 0; i < getNumROI(); i++) {
      kRoi[i].clear();
    }
    kRoi.clear();
  }
}
//
bool G4GMocrenIO::isROIEmpty() {
  if(kRoi.empty()) {
    //if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
    //  G4cout << "!!! ROI data is empty." << G4endl;
    return true;
  } else {
    return false;
  }
}



//----- Track information -----//

int  G4GMocrenIO::getNumTracks() {
  return (int)kSteps.size();
}
int  G4GMocrenIO::getNumTracks4() {
  return (int)kTracks.size();
}
void G4GMocrenIO::addTrack(float * _tracks) {
  kSteps.push_back(_tracks);
}
void G4GMocrenIO::setTracks(std::vector<float *> & _tracks) {
  kSteps = _tracks;
}
std::vector<float *> & G4GMocrenIO::getTracks() {
  return kSteps;
}
void G4GMocrenIO::addTrackColor(unsigned char * _colors) {
  kStepColors.push_back(_colors);
}
void G4GMocrenIO::setTrackColors(std::vector<unsigned char *> & _trackColors) {
  kStepColors = _trackColors;
}
std::vector<unsigned char *> & G4GMocrenIO::getTrackColors() {
  return kStepColors;
}
void G4GMocrenIO::copyTracks(std::vector<float *> & _tracks,
			       std::vector<unsigned char *> & _colors) {
  std::vector<float *>::iterator titr;
  for(titr = kSteps.begin(); titr != kSteps.end(); titr++) {
    float * pts = new float[6];
    for(int i = 0; i < 6; i++) {
      pts[i] = (*titr)[i];
    }
    _tracks.push_back(pts);
  }

  std::vector<unsigned char *>::iterator citr;
  for(citr = kStepColors.begin(); citr != kStepColors.end(); citr++) {
    unsigned char * pts = new unsigned char[3];
    for(int i = 0; i < 3; i++) {
      pts[i] = (*citr)[i];
    }
    _colors.push_back(pts);
  }
}
void G4GMocrenIO::mergeTracks(std::vector<float *> & _tracks,
				std::vector<unsigned char *> & _colors) {
  std::vector<float *>::iterator titr;
  for(titr = _tracks.begin(); titr != _tracks.end(); titr++) {
    addTrack(*titr);
  }

  std::vector<unsigned char *>::iterator citr;
  for(citr = _colors.begin(); citr != _colors.end(); citr++) {
    addTrackColor(*citr);
  }
}
void G4GMocrenIO::addTrack(std::vector<float *> & _steps, unsigned char _color[3]) {

  std::vector<float *>::iterator itr = _steps.begin();
    std::vector<struct GMocrenTrack::Step> steps;
    for(; itr != _steps.end(); itr++) {
      struct GMocrenTrack::Step step;
      for(int i = 0; i < 3; i++) {
	step.startPoint[i] = (*itr)[i];
	step.endPoint[i] = (*itr)[i+3];
      }
      steps.push_back(step);
    }
    GMocrenTrack track;
    track.setTrack(steps);
    track.setColor(_color);
    kTracks.push_back(track);
    
}
void G4GMocrenIO::getTrack(int _num, std::vector<float *> & _steps,
			     std::vector<unsigned char *> & _color) {

  if(_num > (int)kTracks.size()) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "ERROR in getTrack() : " << G4endl;
    G4Exception("G4GMocrenIO::getTrack()",
                "gMocren2003", FatalException,
                "Error.");
  }
  unsigned char * color = new unsigned char[3];
  kTracks[_num].getColor(color);
  _color.push_back(color);

  // steps
  int nsteps = kTracks[_num].getNumberOfSteps();
  for(int isteps = 0; isteps < nsteps; isteps++) {
    float * stepPoints = new float[6];
    kTracks[_num].getStep(stepPoints[0], stepPoints[1], stepPoints[2],
			  stepPoints[3], stepPoints[4], stepPoints[5],
			  isteps);
    _steps.push_back(stepPoints);
  }
}

void G4GMocrenIO::translateTracks(std::vector<float> & _translate) {
  std::vector<class GMocrenTrack>::iterator itr = kTracks.begin();
  for(; itr != kTracks.end(); itr++) {
    itr->translate(_translate);
  }
}




//----- Detector information -----//
int  G4GMocrenIO::getNumberOfDetectors() {
  return (int)kDetectors.size();
}
void G4GMocrenIO::addDetector(std::string & _name,
				std::vector<float *> & _det, 
				unsigned char _color[3]) {

    std::vector<float *>::iterator itr = _det.begin();
    std::vector<struct GMocrenDetector::Edge> edges;
    for(; itr != _det.end(); itr++) {
      struct GMocrenDetector::Edge edge;
      for(int i = 0; i < 3; i++) {
	edge.startPoint[i] = (*itr)[i];
	edge.endPoint[i] = (*itr)[i+3];
      }
      edges.push_back(edge);
    }
    GMocrenDetector detector;
    detector.setDetector(edges);
    detector.setColor(_color);
    detector.setName(_name);
    kDetectors.push_back(detector);
    
}

void G4GMocrenIO::getDetector(int _num, std::vector<float *> & _edges,
				std::vector<unsigned char *> & _color,
				std::string & _detName) {

  if(_num > (int)kDetectors.size()) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "ERROR in getDetector() : " << G4endl;
    
    G4Exception("G4GMocrenIO::getDetector()",
                "gMocren2004", FatalException,
                "Error.");
  }

  _detName = kDetectors[_num].getName();

  unsigned char * color = new unsigned char[3];
  kDetectors[_num].getColor(color);
  _color.push_back(color);

  // edges
  int nedges = kDetectors[_num].getNumberOfEdges();
  for(int ne = 0; ne < nedges; ne++) {
    float * edgePoints = new float[6];
    kDetectors[_num].getEdge(edgePoints[0], edgePoints[1], edgePoints[2],
			     edgePoints[3], edgePoints[4], edgePoints[5],
			     ne);
    _edges.push_back(edgePoints);
  }
}

void G4GMocrenIO::translateDetector(std::vector<float> & _translate) {
  std::vector<class GMocrenDetector>::iterator itr = kDetectors.begin();
  for(; itr != kDetectors.end(); itr++) {
    itr->translate(_translate);
  }
}

// endian conversion
template <typename T>
void G4GMocrenIO::convertEndian(char * _val, T & _rval) {

  if((kLittleEndianOutput && !kLittleEndianInput) ||   // big endian
     (!kLittleEndianOutput && kLittleEndianInput)) {   // little endian

    const int SIZE = sizeof(_rval);
    char ctemp;
    for(int i = 0; i < SIZE/2; i++) {
      ctemp = _val[i];
      _val[i] = _val[SIZE - 1 - i];
      _val[SIZE - 1 - i] = ctemp;
    }
  }
  _rval = *(T *)_val;
}

// inversion of byte order
template <typename T>
void G4GMocrenIO::invertByteOrder(char * _val, T & _rval) {

  const int SIZE = sizeof(_rval);
  //char * cval = new char[SIZE];
  union {
    char cu[16];
    T tu;
  } uni;
  for(int i = 0; i < SIZE; i++) {
    uni.cu[i] = _val[SIZE-1-i];
    //cval[i] = _val[SIZE-i-1];
  }
  //_rval = *(T *)cval;
  _rval = uni.tu;
  //delete [] cval;
}

//----- kVerbose information -----//
void G4GMocrenIO::setVerboseLevel(int _level) {
  kVerbose = _level;
}

