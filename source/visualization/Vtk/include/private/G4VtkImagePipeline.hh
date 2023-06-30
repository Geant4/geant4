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

#ifndef G4VTKIMAGEPIPELINE_HH
#define G4VTKIMAGEPIPELINE_HH

#include "G4Transform3D.hh"
#include "G4Types.hh"
#include "G4VVtkPipeline.hh"

#include "vtkImageData.h"
#include "vtkImageReader2.h"
#include "vtkImageReader2Factory.h"
#include "vtkSmartPointer.h"

#include <vector>

class vtkImageAlgorithm;
class vtkImageActor;

class G4String;
class G4VtkVisContext;

class G4VtkImagePipeline : G4VVtkPipeline
{
  public:
    G4VtkImagePipeline(const G4String& name, const G4VtkVisContext& vc);
    ~G4VtkImagePipeline() override = default;
    void AddFilter(vtkSmartPointer<vtkImageAlgorithm> f) { filters.push_back(f); }
    vtkSmartPointer<vtkImageAlgorithm> GetFilter(G4int iFilter) { return filters[iFilter]; }
    vtkSmartPointer<vtkImageAlgorithm> GetFinalFilter() { return *filters.end(); }

    void AddChildPipeline(G4VtkImagePipeline* child)
    {
      childPipelines.push_back(child);
      if (child->GetDisableParent()) {
        Disable();
      }
    }
    G4VtkImagePipeline* GetChildPipeline(G4int iPipeline) { return childPipelines[iPipeline]; }

    virtual void SetImage(const G4String& fileName);
    virtual void SetTransformation(const G4Transform3D& transformation);

    void Print() override;
    void Modified() override;
    void Clear() override;

    virtual vtkSmartPointer<vtkImageActor> GetActor() { return actor; }

    void Disable() override;
    void Enable() override;

  private:
    std::vector<vtkSmartPointer<vtkImageAlgorithm>>
      filters;  // derived types can store filters in this vector
    std::vector<G4VtkImagePipeline*> childPipelines;
    vtkSmartPointer<vtkImageActor> actor;  // all pipelines require an actor
};

#endif  // G4VTKIMAGEPIPELINE_HH
