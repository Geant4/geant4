
#include <vector>

#include "interfaces/IPlotter.h"

#include "interfaces/IHistogram.h"
#include "interfaces/IHistogram1D.h"
#include "interfaces/IHistogram2D.h"

class IHistoFactory;
class IVectorFactory;
class IPlotter;

class XrayTelHistogram {

public:
  XrayTelHistogram();
  virtual ~XrayTelHistogram();

  bool book();
  bool finish();

  bool plot(IHistogram1D *);
  bool plot(IHistogram2D *);

  IPlotter * getPlotter() { return pl; }

public:
  // to have access to the histograms 
  vector<IHistogram1D *> * getH1DList() { return &h1dList; }
  vector<IHistogram2D *> * getH2DList() { return &h2dList; }

private:
  IHistoFactory  *hFact;
  IVectorFactory *vFact;
  IPlotter       *pl;

  vector<IHistogram1D *> h1dList;
  vector<IHistogram2D *> h2dList;

};
