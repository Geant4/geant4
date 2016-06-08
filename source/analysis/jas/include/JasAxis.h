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
#ifndef JASAXIS_H
#define JASAXIS_H

#if defined(G4ANALYSIS_BUILD_JAS)

#include <IAxis.h>

class JasAxis : public IAxis {
public:
  virtual ~JasAxis() {}
public:
  virtual double lowerEdge() const;
  virtual double upperEdge() const;
  virtual int bins() const;
  virtual double binLowerEdge(int) const;
  virtual double binUpperEdge(int) const;
  virtual double binWidth(int) const;
  virtual double binCentre(int) const;
  virtual int coordToIndex(double) const;
};

#endif

#endif
