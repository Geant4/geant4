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
#ifndef JASHISTOGRAM2D_H
#define JASHISTOGRAM2D_H

#if defined(G4ANALYSIS_BUILD_JAS)

#include <IHistogram2D.h>

#include "JasAxis.h"

class JasHistogramFactory;

class JasHistogram2D : public IHistogram2D {
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
public: // IHistogram2D IMethods :
  virtual void fill( double x, double y, double weight = 1 );
  virtual int binEntries( int indexX, int indexY ) const;
  virtual int binEntriesX( int indexX ) const;
  virtual int binEntriesY( int indexY ) const;
  virtual double binHeight( int indexX, int indexY ) const;
  virtual double binHeightX( int indexX ) const;
  virtual double binHeightY( int indexY ) const;
  virtual double binError( int indexX, int indexY ) const;
  virtual double meanX() const;
  virtual double meanY() const;
  virtual double rmsX() const;
  virtual double rmsY() const;
  virtual int minBinX() const;
  virtual int minBinY() const;
  virtual int maxBinX() const;
  virtual int maxBinY() const;
  virtual const IAxis& xAxis() const;
  virtual const IAxis& yAxis() const;
  virtual int coordToIndexX( double coordX ) const;
  virtual int coordToIndexY( double coordY ) const;
  virtual IHistogram1D* projectionX() const;
  virtual IHistogram1D* projectionY() const;
  virtual IHistogram1D* sliceX( int indexY ) const;
  virtual IHistogram1D* sliceY( int indexX ) const;
  virtual IHistogram1D* sliceX( int indexY1, int indexY2 ) const;
  virtual IHistogram1D* sliceY( int indexX1, int indexX2 ) const;
public:
 JasHistogram2D(JasHistogramFactory *factory, const char *title);
 virtual ~JasHistogram2D();
private:
 jmethodID fillMethod;
 jobject fJasHist;
 JNIEnv* fEnv;
 JasAxis fAxis;  
 std_string fTitle;
};

#endif

#endif
