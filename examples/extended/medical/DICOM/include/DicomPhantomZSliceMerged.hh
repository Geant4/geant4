//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DicomPhantomZSliceMerged.hh 76941 2013-11-19 09:55:03Z gcosmo $
//
/// \file medical/DICOM/include/DicomPhantomZSliceMerged.hh
/// \brief Definition of the DicomPhantomZSliceMerged class
//
//
// The code was written by :
//      * Jonathan Madsen : jonathan.madsen@cern.ch (12/18/2012)
//
//
//*******************************************************


#ifndef dicomphantomzslicemerged_hh_
#define dicomphantomzslicemerged_hh_

#include <map>
#include "globals.hh"

#include "DicomPhantomZSliceHeader.hh"

class DicomPhantomZSliceMerged
{
public:
    // Constructor and Destructors
    DicomPhantomZSliceMerged();
    ~DicomPhantomZSliceMerged();

public:
    // Public functions
    void AddZSlice(DicomPhantomZSliceHeader* val) {
        fSlices[val->GetSliceLocation()] = val;
        //std::cout << "Slice Location : " << val->GetSliceLocation()/mm << " mm" << std::endl;
    }

    void CheckSlices();

    inline void DumpExcessMemory();

private:
    // Private functions

private:
    // Private variables
    std::map<G4double,DicomPhantomZSliceHeader*> fSlices;


};

inline void DicomPhantomZSliceMerged::DumpExcessMemory()
{
    for(std::map<G4double,DicomPhantomZSliceHeader*>::iterator ite = fSlices.begin();
        ite != fSlices.end(); ++ite) {
        ite->second->DumpExcessMemory();
    }
}


#endif
