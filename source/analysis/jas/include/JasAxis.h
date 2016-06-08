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
