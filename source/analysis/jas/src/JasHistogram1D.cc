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
// JasHistogram.cpp: implementation of the JasHistogram1D class.
//////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_BUILD_JAS

#include <jni.h>
#include <string.h>
#include <stdlib.h>
#include "JasHistogram1D.h"
#include "JasHistogramFactory.h"

JasHistogram1D::JasHistogram1D(JasHistogramFactory* factory,  const char* title)
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

JasHistogram1D::~JasHistogram1D()
{
  if(fEnv==NULL) return;
  fEnv->DeleteLocalRef(fJasHist);
}

void JasHistogram1D::fill(double value,double)
{
  if(fEnv==NULL) return;
  fEnv->CallVoidMethod(fJasHist,fillMethod,value);
}

// Today dummy methods....
#define IWeight     double
#define ICoordinate double
#define ISize       int
#define IBin        int

// IHistogram :
std_string JasHistogram1D::title() const {return fTitle;}
int JasHistogram1D::dimensions() const {return 1;}
ISize JasHistogram1D::entries() const {return 0;}
ISize JasHistogram1D::allEntries() const {return 0;}
ISize JasHistogram1D::extraEntries() const {return 0;}
IWeight JasHistogram1D::equivalentBinEntries() const {return 0;}
IWeight JasHistogram1D::sumBinHeights() const {return 0;}
IWeight JasHistogram1D::sumAllBinHeights() const {return 0;}
IWeight JasHistogram1D::sumExtraBinHeights() const {return 0;}
IWeight JasHistogram1D::minBinHeight() const {return 0;}
IWeight JasHistogram1D::maxBinHeight() const {return 0;}
void JasHistogram1D::reset() {}
IAnnotation* JasHistogram1D::annotation() {return 0;}
// IHistogram1D :
IBin JasHistogram1D::minBin() const {return 0;}
IBin JasHistogram1D::maxBin() const {return 0;}
IBin JasHistogram1D::coordToIndex(ICoordinate) const {return 0;}
ICoordinate JasHistogram1D::mean() const {return 0;}
ICoordinate JasHistogram1D::rms() const {return 0;}

// Bins :
ISize JasHistogram1D::binEntries(IBin) const {return 0;}
IWeight JasHistogram1D::binHeight(IBin) const {return 0;}
IWeight JasHistogram1D::binError(IBin) const {return 0;}

const IAxis& JasHistogram1D::xAxis() const {return fAxis;}
  
#endif
