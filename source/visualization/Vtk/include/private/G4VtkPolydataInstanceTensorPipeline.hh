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

#ifndef G4VTKPOLYDATAINSTANCETENSORPIPELINE_HH
#define G4VTKPOLYDATAINSTANCETENSORPIPELINE_HH

#include "G4VtkPolydataInstancePipeline.hh"

class G4String;
class G4VtkVisContext;
class vtkTensorGlyphColor;

class G4VtkPolydataInstanceTensorPipeline : public G4VtkPolydataInstancePipeline
{
  public:
    G4VtkPolydataInstanceTensorPipeline(G4String name, const G4VtkVisContext& vc);
    ~G4VtkPolydataInstanceTensorPipeline() override = default;

    void Print() override;


    static std::size_t MakeHash(const G4Polyhedron& p, const G4VtkVisContext& vc);

  protected:
    vtkSmartPointer<vtkPolyData> instancePolydata;
    vtkSmartPointer<vtkTensorGlyphColor> instanceTensorGlyph;
};

#endif  // G4VTKPOLYDATAINSTANCETENSORPIPELINE_HH
