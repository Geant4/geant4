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
//
// $Id: DataHists.hh,v 1.2 2003-05-27 14:35:04 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef DataHists_h
#define DataHists_h 1

#ifdef USE_AIDA
#include "AIDA/AIDA.h"
#include "globals.hh"

using namespace AIDA;

class EnergyAngleCrossSection;


class DataHists
{
  public:
    DataHists();
    DataHists(G4double lo, G4double hi, G4double size);
    ~DataHists();

    void CreateHists(IHistogramFactory*);
    void PlotData(EnergyAngleCrossSection*);
    
  private:

    IHistogramFactory* hFactory;

    // pi energy distribution histograms

    G4std::vector<IHistogram1D*> datahists;

    G4double loBinEdge;
    G4double hiBinEdge;
    G4double binSize;

};
#endif

#endif

    




