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
// JasHistogram.cpp: implementation of the JasHistogram2D class.
//////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_BUILD_JAS

#include <jni.h>
#include <string.h>
#include <stdlib.h>

#include "JasHistogram2D.h"
#include "JasHistogramFactory.h"

JasHistogram2D::JasHistogram2D(JasHistogramFactory* factory,  const char* title)
:fillMethod(NULL)
,fJasHist(NULL)
,fEnv(NULL)
{
  fEnv =  factory->getEnv();
  if(fEnv==NULL) return;
  // Create the corresponding Java histogram
  jclass cls = fEnv->FindClass("hep/analysis/Histogram");
  if (cls == NULL) {
    factory->error("Could not create hep.analysis.Histogram");
    fEnv = 0;
    return;
  }
  jmethodID constructor = 
    fEnv->GetMethodID(cls, "<init>", "(Ljava/lang/String;)V");
  if (constructor == NULL) {
    factory->error("Could not find constructor for hep.analysis.Histogram");
    fEnv = 0;
    return;
  }
  fillMethod = fEnv->GetMethodID(cls,"fill","(D)V");
  if (fillMethod == NULL) {
    factory->error("Could not find fill method for hep.analysis.Histogram");
    fEnv = 0;
    return;
  }
  jstring jtitle = fEnv->NewStringUTF(title);
  fJasHist = fEnv->NewObject(cls,constructor,jtitle);
}

JasHistogram2D::~JasHistogram2D()
{
  if(fEnv==NULL) return;
  fEnv->DeleteLocalRef(fJasHist);
}

void JasHistogram2D::fill(double x,double y, double w)
{
  if(fEnv==NULL) return;
  fEnv->CallVoidMethod(fJasHist,fillMethod,x);
}

// IHistogram :
std_string JasHistogram2D::title() const {return fTitle;}
int JasHistogram2D::dimensions() const {return 1;}
int JasHistogram2D::entries() const {return 0;}
int JasHistogram2D::allEntries() const {return 0;}
int JasHistogram2D::extraEntries() const {return 0;}
double JasHistogram2D::equivalentBinEntries() const {return 0;}
double JasHistogram2D::sumBinHeights() const {return 0;}
double JasHistogram2D::sumAllBinHeights() const {return 0;}
double JasHistogram2D::sumExtraBinHeights() const {return 0;}
double JasHistogram2D::minBinHeight() const {return 0;}
double JasHistogram2D::maxBinHeight() const {return 0;}
void JasHistogram2D::reset() {}
IAnnotation* JasHistogram2D::annotation() {return 0;}

// IHistogram2D :
int JasHistogram2D::binEntries( int indexX, int indexY ) const {return 0;}
int JasHistogram2D::binEntriesX( int indexX ) const {return 0;}
int JasHistogram2D::binEntriesY( int indexY ) const {return 0;}
double JasHistogram2D::binHeight( int indexX, int indexY ) const {return 0;}
double JasHistogram2D::binHeightX( int indexX ) const {return 0;}
double JasHistogram2D::binHeightY( int indexY ) const {return 0;}
double JasHistogram2D::binError( int indexX, int indexY ) const {return 0;}
double JasHistogram2D::meanX() const {return 0;}
double JasHistogram2D::meanY() const {return 0;}
double JasHistogram2D::rmsX() const {return 0;}
double JasHistogram2D::rmsY() const {return 0;}
int JasHistogram2D::minBinX() const {return 0;}
int JasHistogram2D::minBinY() const {return 0;}
int JasHistogram2D::maxBinX() const {return 0;}
int JasHistogram2D::maxBinY() const {return 0;}
const IAxis& JasHistogram2D::xAxis() const {return fAxis;}
const IAxis& JasHistogram2D::yAxis() const {return fAxis;}
int JasHistogram2D::coordToIndexX( double coordX ) const {return 0;}
int JasHistogram2D::coordToIndexY( double coordY ) const {return 0;}
IHistogram1D* JasHistogram2D::projectionX() const {return 0;}
IHistogram1D* JasHistogram2D::projectionY() const {return 0;}
IHistogram1D* JasHistogram2D::sliceX( int indexY ) const {return 0;}
IHistogram1D* JasHistogram2D::sliceY( int indexX ) const {return 0;}
IHistogram1D* JasHistogram2D::sliceX( int indexY1, int indexY2 ) const {return 0;}
IHistogram1D* JasHistogram2D::sliceY( int indexX1, int indexX2 ) const {return 0;}

#endif
