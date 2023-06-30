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
// File I/O manager class for writing or reading calcuated dose
// distribution and some event information
//
//
//  Mar. 31, 2009 :  release for the gMocrenFile driver
//
//                               Akinori Kimura
//
//                               gMocren home page:
//                               http://geant4.kek.jp/gMocren/
//
#ifndef GMOCRENIO_HH
#define GMOCRENIO_HH

#include <vector>
#include <string>
#include <fstream>
#include <map>

//
//----- GMocrenDataPrimitive class -----//
// data primitive class for volume data
//
template <typename T> class GMocrenDataPrimitive {
protected:
  int kSize[3];
  double kScale;
  T kMinmax[2];
  float kCenter[3];
  std::vector<T *> kImage;
  std::string kDataName;
  //std::vector<std::vector<T>> image;

public:
  GMocrenDataPrimitive();
  //GMocrenDataPrimitive(GMocrenDataPrimitive<T> & _prim);
  ~GMocrenDataPrimitive();

  GMocrenDataPrimitive<T> & operator = (const GMocrenDataPrimitive<T> & _right);
  GMocrenDataPrimitive<T> & operator + (const GMocrenDataPrimitive<T> & _right);
  GMocrenDataPrimitive<T> & operator += (const GMocrenDataPrimitive<T> & _right);

  void clear();
  void clearImage();
  void setSize(int _size[3]);
  void getSize(int _size[3]);
  void setScale(double & _scale);
  double getScale();
  void setMinMax(T _minmax[2]);
  void getMinMax(T _minmax[2]);
  void setImage(std::vector<T *> & _image);
  void addImage(T * _image);
  std::vector<T *> & getImage();
  T * getImage(int _z);  // get image of each layer
  void setCenterPosition(float _center[3]);
  void getCenterPosition(float _center[3]);
  void setName(std::string & _name);
  std::string getName();
};


//
//----- GMocrenTrack class -----//
//
class GMocrenTrack {
public:
  struct Step {
    float startPoint[3];
    float endPoint[3];
  };
protected:
  std::vector<struct Step> kTrack;
  unsigned char kColor[3];

public:
  GMocrenTrack();
  ~GMocrenTrack(){;}

  
  int getNumberOfSteps() {return (int)kTrack.size();}
  void addStep(float _startx, float _starty, float _startz,
	       float _endx, float _endy, float _endz);
  void getStep(float & _startx, float & _starty, float & _startz,
	       float & _endx, float & _endy, float & _endz,
	       int _num);
  void setTrack(std::vector<struct Step> & _aTrack) {kTrack = _aTrack;}
  void setColor(unsigned char _color[3]) {
    for(int i = 0; i < 3; i++) kColor[i] = _color[i];
  }
  void getColor(unsigned char _color[3]) {
    for(int i = 0; i < 3; i++) _color[i] = kColor[i];
  }
  void translate(std::vector<float> & _tranlate);
};



//
//----- GMocrenDetector class -----//
//
class GMocrenDetector {
public:
  struct Edge {
    float startPoint[3];
    float endPoint[3];
  };
protected:
  std::vector<struct Edge> kDetector;
  unsigned char kColor[3];
  std::string kName;

public:
  GMocrenDetector();
  ~GMocrenDetector(){;}

  
  int getNumberOfEdges() {return (int)kDetector.size();}
  void addEdge(float _startx, float _starty, float _startz,
	       float _endx, float _endy, float _endz);
  void getEdge(float & _startx, float & _starty, float & _startz,
	       float & _endx, float & _endy, float & _endz,
	       int _num);
  void setDetector(std::vector<struct Edge> & _aDetector) {kDetector = _aDetector;}
  void setColor(unsigned char _color[3]) {
    for(int i = 0; i < 3; i++) kColor[i] = _color[i];
  }
  void getColor(unsigned char _color[3]) {
    for(int i = 0; i < 3; i++) _color[i] = kColor[i];
  }
  void setName(std::string & _name) { kName = _name;}
  std::string getName() {return kName;}

  void translate(std::vector<float> & _tranlate);
};


//
//----- G4GMocrenIO class -----//
//
class G4GMocrenIO {
public:
  // file id
  static std::string kId;

  // file version
  static std::string kVersion;

  // data file name 
  static std::string kFileName;

  // file data endian: little or not
  static char kLittleEndianInput;
  static char kLittleEndianOutput;

  static std::string kComment;

  // number of events
  static int kNumberOfEvents;

  // pointer to the modality image data
  static unsigned int kPointerToModalityData;
  // pointer to the dose distribution image data
  static std::vector<unsigned int> kPointerToDoseDistData;
  // pointer to the ROI image data
  static unsigned int kPointerToROIData;
  // pointer to the track data
  static unsigned int kPointerToTrackData;
  // pointer to the detector data
  static unsigned int kPointerToDetectorData;

  // voxel spacing (universal size)
  static float kVoxelSpacing[3];

  //----- modality image -----//
  static class GMocrenDataPrimitive<short> kModality;
  // density map to modality (CT) values
  static std::vector<float> kModalityImageDensityMap;
  static std::string kModalityUnit;

  //----- dose distribution -----//
  static std::vector<class GMocrenDataPrimitive<double> > kDose;
  //std::vector<short *> kShortDose;
  static std::string kDoseUnit;

  //----- RoI -----//
  static std::vector<class GMocrenDataPrimitive<short> > kRoi;

  //----- track information -----//
  static std::vector<float *> kSteps; // begin (x,y,z), end (x,y,z)
  static std::vector<unsigned char *> kStepColors; // r, g, b

  static std::vector<class GMocrenTrack> kTracks;
  bool kTracksWillBeStored;

  //----- detector information -----//
  static std::vector<class GMocrenDetector> kDetectors;

  //----- verbose information -----//
  static int kVerbose; // verbose level :  0 - 5 (none - overtalk)

public:
  // constructor
  G4GMocrenIO();
  // destructor
  ~G4GMocrenIO();

  // initialize
  void initialize();

  // set the gMocren data file name
  void setFileName(std::string & _filename) {kFileName = _filename;}
  void setFileName(char * _filename) {kFileName = _filename;}
  // get the gMocren data file name
  std::string & getFileName() {return kFileName;}
  // store all data in the gMocren data file
  bool storeData(char * _filename); // interface for version 4
  bool storeData();
  bool storeData2(char * _filename); // version 2
  bool storeData2();
  bool storeData3(char * _filename); // version 3
  bool storeData3();
  bool storeData4(char * _filename); // version 4
  bool storeData4();
  // retrieve all data from the gMocren data file
  bool retrieveData(char * _filename); // interface
  bool retrieveData();
  bool retrieveData2(char * _filename); //version 2
  bool retrieveData2();
  bool retrieveData3(char * _filename); // version 3
  bool retrieveData3();
  bool retrieveData4(char * _filename); // version 4
  bool retrieveData4();
    
  // get & set the file id
  std::string & getID() {return kId;}
  void setID();
  void setID(std::string & _id) {kId = _id;}

  // get & set the file version
  std::string & getVersion();
  void setVersion(std::string & _version);

  // set endians of input/output data
  void setLittleEndianInput(bool _little);
  void setLittleEndianOutput(bool _little);

  // get & set comment
  std::string & getComment() {return kComment;}
  void setComment(std::string & _comment) {kComment = _comment;}
  

  // voxel spacing
  void setVoxelSpacing(float _spacing[3]);
  void getVoxelSpacing(float _spacing[3]);

  // get & set number of events
  int & getNumberOfEvents();
  void setNumberOfEvents(int & _numberOfEvents);
  void addOneEvent();

  // set pointer the modality image data
  void setPointerToModalityData(unsigned int & _pointer);
  unsigned int getPointerToModalityData();
  // set pointer the dose distribution image data
  void addPointerToDoseDistData(unsigned int & _pointer);
  unsigned int getPointerToDoseDistData(int _elem = 0);
  // set pointer the ROI image data
  void setPointerToROIData(unsigned int & _pointer);
  unsigned int getPointerToROIData();
  // set pointer the track data
  void setPointerToTrackData(unsigned int & _pointer);
  unsigned int getPointerToTrackData();
private:
  // calculate pointers
  void calcPointers4();
  void calcPointers3();
  void calcPointers2();


  //----- Modality image -----//
public:
  // get & set the modality image size
  void getModalityImageSize(int _size[3]);
  void setModalityImageSize(int _size[3]);
  // get & set the modality image spacing size
  void getModalityImageVoxelSpacing(float _size[3]); // un-usable
  void setModalityImageVoxelSpacing(float _size[3]); // un-usable
  // get & set the modality image size
  void setModalityImageScale(double & _scale);
  double getModalityImageScale();
  // set the modality image in CT 
  void setModalityImage(short * _image);
  short * getModalityImage(int _z);
  void clearModalityImage();
  // set/get the modality image density map
  void setModalityImageDensityMap(std::vector<float> & _map);
  std::vector<float> & getModalityImageDensityMap();
  // set the modality image min./max.
  void setModalityImageMinMax(short _minmax[2]);
  // get min. & max. of the modality image 
  void getModalityImageMinMax(short _minmax[2]);
  short getModalityImageMax();
  short getModalityImageMin();
  // set center of the modality image position
  void setModalityCenterPosition(float _center[3]);
  void getModalityCenterPosition(float _center[3]);
  // get & set the modality image unit
  std::string getModalityImageUnit();
  void setModalityImageUnit(std::string & _unit);

  short convertDensityToHU(float & _dens);

  //----- Dose distribution -----//

  // instanciate a dose distribution data object
  void newDoseDist();
  // get number of dose distribion data
  int getNumDoseDist();
  // get & set the dose distribution unit
  std::string getDoseDistUnit(int _num = 0);
  void setDoseDistUnit(std::string & _unit, int _num = 0);
  // get & set the dose distribution image size
  void getDoseDistSize(int _size[3], int _num = 0);
  void setDoseDistSize(int _size[3], int _num = 0);
  // get min. & max. of the dose distribution image
  void setDoseDistMinMax(short _minmax[2], int _num = 0);
  void getDoseDistMinMax(short _minmax[2], int _num = 0);
  // get min. & max. of the dose distribution 
  void setDoseDistMinMax(double _minmax[2], int _num = 0);
  void getDoseDistMinMax(double _minmax[2], int _num = 0);
  // get & set scale value of the dose distribution for the image
  void setDoseDistScale(double & _scale, int _num = 0);
  double getDoseDistScale(int _num = 0);
  // set the dose distribution image
  void setShortDoseDist(short * _image, int _num = 0);
  void getShortDoseDist(short * _data, int _z, int _num = 0);
  void getShortDoseDistMinMax(short _minmax[2], int _num = 0);
  // set the dose distribution 
  void setDoseDist(double * _image, int _num = 0);
  double * getDoseDist(int _z, int _num = 0);
  // add another dose ditribution map to this map
  bool addDoseDist(std::vector<double *> & _image, int _num = 0);

  // get & get center position of calculated dose region
  void getDoseDistCenterPosition(float _center[3], int _num = 0);
  void setDoseDistCenterPosition(float _center[3], int _num = 0);

  // get & get name of calculated dose distribution
  std::string getDoseDistName(int _num = 0);
  void setDoseDistName(std::string _name, int _num = 0);

  // copy dose distributions
  void copyDoseDist(std::vector<class GMocrenDataPrimitive<double> > & _dose);
  // merge two dose distributions
  bool mergeDoseDist(std::vector<class GMocrenDataPrimitive<double> > & _dose);

  // clear all dose distributions
  void clearDoseDistAll();
protected:
  // check whether dose variable is empty or not
  bool isDoseEmpty();
  // calcuated scale value to convert dose distribution into image
  void calcDoseDistScale();

public:
  //----- RoI -----//

  // instanciate an RoI data object
  void newROI();
  // get number of RoI data
  int getNumROI();
  // get & set the ROI image scale
  double getROIScale(int _num = 0);
  void setROIScale(double & _scale, int _num = 0);
  // get & set the ROI image 
  short * getROI(int _z, int _num = 0);
  void setROI(short * _image, int _num = 0);
  // get & set the ROI image size
  void getROISize(int _size[3], int _num = 0);
  void setROISize(int _size[3], int _num = 0);
  // get & set position of the ROI region center
  void getROICenterPosition(float _center[3], int _num = 0);
  void setROICenterPosition(float _center[3], int _num = 0);
  // get & set the ROI image min. and max.
  void getROIMinMax(short _minmax[2], int _num = 0);
  void setROIMinMax(short _minmax[2], int _num = 0);
  void clearROIAll();
protected:
  // check whether RoI variable is empty or not
  bool isROIEmpty();


public:
  //----- Track -----//
  // get number of tracks
  int getNumTracks();
  int getNumTracks4();
  // get & set tracks
  std::vector<float *> & getTracks();
  void getTrack(int _num, std::vector<float *> & _steps, 
		std::vector<unsigned char * > & _color);
  void addTrack(float * _tracks);
  void setTracks(std::vector<float *> & _tracks);
  std::vector<unsigned char *> & getTrackColors();
  void addTrackColor(unsigned char * _colors);
  void setTrackColors(std::vector<unsigned char *> & _trackColors);
  void copyTracks(std::vector<float *> & _tracks, std::vector<unsigned char *> & _colors);
  void mergeTracks(std::vector<float *> & _tracks, std::vector<unsigned char *> & _colors);
  void addTrack(std::vector<float *> & _steps, unsigned char _color[3]);

  void notStoredTracks() {kTracksWillBeStored = false;};
  void translateTracks(std::vector<float> & _translateo);
  void clearTracks() {kTracks.clear();}


  //----- Detectors -----//
  // get number of detectors
  int getNumberOfDetectors();
  // add one detector which consists of edges (float[6])
  void addDetector(std::string & _name, std::vector<float *> & _det, unsigned char _color[3]);
  void getDetector(int _num, std::vector<float *> & _edges,
		   std::vector<unsigned char *> & _color,
		   std::string & _detectorName);
  void translateDetector(std::vector<float> & _translate);
  void clearDetector() {kDetectors.clear();}

protected:
  // endian conversion
  template <typename Type> void convertEndian(char *, Type &);
  // byte order inversion
  template <typename T> void invertByteOrder(char * _val, T & _rval);


public:
  //----- verbose information -----//
  void setVerboseLevel(int _level);

};

#endif

