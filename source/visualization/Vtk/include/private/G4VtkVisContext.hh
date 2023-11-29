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

#ifndef G4VTKVISCONTEXT_HH
#define G4VTKVISCONTEXT_HH

#include "G4ViewParameters.hh"

class G4VtkViewer;
class G4VisAttributes;

class G4VtkVisContext
{
  public:
    const G4VtkViewer* fViewer;
    const G4VisAttributes* fVisAtt;
    G4bool fProcessing2D;
    const G4Transform3D& fTransform;
    G4int fDepth = -1;
    G4String fDescription;
    G4double red, green, blue, alpha;

    G4double fSize;
    G4ViewParameters::DrawingStyle fDrawingStyle;

    G4VtkVisContext()
      : fViewer(nullptr),
        fVisAtt(nullptr),
        fProcessing2D(false),
        fTransform(G4Transform3D::Identity),
        fDepth(0),
        fDescription(""),
        red(0),
        green(0),
        blue(0),
        alpha(0)
    {}

    G4VtkVisContext(const G4VtkViewer* viewer, const G4VisAttributes* visAtt, bool processing2D,
                    const G4Transform3D& transform)
      : fViewer(viewer), fVisAtt(visAtt), fProcessing2D(processing2D), fTransform(transform)
    {}

    G4VtkVisContext(const G4VtkVisContext& vc2)
      : fViewer(vc2.fViewer),
        fVisAtt(vc2.fVisAtt),
        fProcessing2D(vc2.fProcessing2D),
        fTransform(vc2.fTransform),
        fDepth(vc2.fDepth),
        fDescription(vc2.fDescription),
        red(vc2.red),
        green(vc2.green),
        blue(vc2.blue),
        alpha(vc2.alpha),
        fDrawingStyle(vc2.fDrawingStyle)
    {}
};

#endif
