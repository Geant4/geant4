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
#ifndef JASHISTOGRAM1D_H
#define JASHISTOGRAM1D_H

#if defined(G4ANALYSIS_BUILD_JAS)

#include <IHistogram1D.h>

#include "JasAxis.h"

class JasHistogramFactory;

class JasHistogram1D : public IHistogram1D {
public: // IHistogram IMethods :
  virtual std_string title() const;
  virtual int dimensions() const;
  virtual int entries() const;
  virtual int allEntries() const; 
  virtual int extraEntries() const;
  virtual double equivalentBinEntries() const;
  virtual double sumBinHeights() const;
  virtual double sumAllBinHeights() const;
  virtual double sumExtraBinHeights() const;
  virtual double minBinHeight() const;
  virtual double maxBinHeight() const;
  virtual void reset();
  virtual IAnnotation* annotation();
public: // IHistogram1D IMethods :
  // The fill :
  virtual void fill(double,double = 1.0);

  // Partition :
  virtual int minBin() const;
  virtual int maxBin() const;
  virtual int coordToIndex(double) const;

  virtual double mean() const;
  virtual double rms() const;

  // Bins :
  virtual int binEntries(int) const;
  virtual double binHeight(int) const;
  virtual double binError(int) const;

  // Axis :  
  virtual const IAxis& xAxis() const;
public:
 JasHistogram1D(JasHistogramFactory *factory, const char *title);
 virtual ~JasHistogram1D();
private:
 jmethodID fillMethod;
 jobject fJasHist;
 JNIEnv* fEnv;
 JasAxis fAxis;  
 std_string fTitle;
};

#endif

#endif
