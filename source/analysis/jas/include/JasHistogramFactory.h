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
// JHistogramFactory.h: interface for the JHistogramFactory class.
//////////////////////////////////////////////////////////////////////
#ifndef JASHISTOGRAMFACTORY
#define JASHISTOGRAMFACTORY

#if defined(G4ANALYSIS_BUILD_JAS)

#include <jni.h>

#include <IHistogramFactory.h>

class JasHistogramFactory  
: public IHistogramFactory
{
 public: // IMethods :
  virtual IHistogram1D* create1D(const std_string&,
				 int=0,double=0,double=0,
				 int=0);
  virtual IHistogram2D* create2D(const std_string&,
				 int=0,double=0,double=0,
				 int=0,double=0,double=0,
				 int=0);
  virtual void destroy(IHistogram*);
 public:
  JasHistogramFactory();
  virtual ~JasHistogramFactory();
  JNIEnv* getEnv();
  void error(const char* message);
 private:
  jobject fJob;
  void createJVM();
  void destroyJVM();
};

#endif

#endif
