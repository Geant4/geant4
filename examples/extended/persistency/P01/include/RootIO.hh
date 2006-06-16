// $Id: RootIO.hh,v 1.1 2006-06-16 12:31:12 gcosmo Exp $
#ifndef INCLUDE_ROOTIO_HH 
#define INCLUDE_ROOTIO_HH 1

// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"

#include "ExP01TrackerHit.hh"

/** @class rootio rootio.hh include/rootio.hh
 *   
 *
 *  @author Witold POKORSKI
 *  @date   2005-10-27
 */

class RootIO 
{
public: 
  virtual ~RootIO();
  
  static RootIO* GetInstance();
  void Write(std::vector<ExP01TrackerHit*>*);

protected:
  RootIO(); 
  
private:

  TFile* fo;
  int Nevents;
  
};
#endif // INCLUDE_ROOTIO_HH
