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
//
/// \file medical/DICOM/include/DicomPhantomParameterisationColour.hh
/// \brief Definition of the DicomPhantomParameterisationColour class
//

#ifndef DicomPhantomParameterisationColour_HH
#define DicomPhantomParameterisationColour_HH

#include <map>

#include "G4PhantomParameterisation.hh"
class G4VisAttributes;

// *********************************************************************
/// Class inherited from G4PhantomParameterisation to provide different
//  colour for each material
///
/// History:
/// - Created.    5 December 2007
/// \author       P. Arce
// *********************************************************************

class DicomPhantomParameterisationColour : public G4PhantomParameterisation
{
public:
    typedef std::map<G4String,G4VisAttributes*> ColourMap_t;

    static G4String defaultColorFile;

public:  // with description
    DicomPhantomParameterisationColour(G4String colorFile =
                                       defaultColorFile);
    ~DicomPhantomParameterisationColour();

    virtual G4Material* ComputeMaterial(const G4int repNo,
                                        G4VPhysicalVolume *currentVol,
                                        const G4VTouchable *parentTouch=0);

    const ColourMap_t& GetColourMap() const { return fColours; }
    ColourMap_t& GetColourMap() { return fColours; }

private:
    void ReadColourData(G4String colourFile);

private:
    ColourMap_t fColours;
    std::map<G4int, G4VisAttributes*> mColours;
};


#endif
