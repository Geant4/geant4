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
// $Id: RootIO.hh,v 1.2 2006-06-28 13:58:01 gunter Exp $
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
