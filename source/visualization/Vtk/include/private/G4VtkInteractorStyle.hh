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
// Created by Stewart Boogert on 28/02/2023.
//

#ifndef G4VTKINTERACTORSTYLE_HH
#define G4VTKINTERACTORSTYLE_HH

#include "G4ios.hh"
#include "vtkInteractionStyleModule.h"  // For export macro
#include "vtkInteractorStyleMultiTouchCamera.h"
#include <vtkObjectFactory.h>

// Define interaction style
class VTKINTERACTIONSTYLE_EXPORT G4VtkInteractorStyle : public vtkInteractorStyleMultiTouchCamera
{
  public:
    static G4VtkInteractorStyle* New();
    vtkTypeMacro(G4VtkInteractorStyle, vtkInteractorStyleMultiTouchCamera)

    void OnMouseMove() override
    {
      vtkInteractorStyleMultiTouchCamera::OnMouseMove();
    }

    void OnLeftButtonDown() override
    {
      // Forward events
      vtkInteractorStyleMultiTouchCamera::OnLeftButtonDown();
    }


    void OnLeftButtonUp() override
    {
      // Forward events
      vtkInteractorStyleMultiTouchCamera::OnLeftButtonUp();
    }

    void OnMiddleButtonDown() override
    {
      // Forward events
      vtkInteractorStyleMultiTouchCamera::OnMiddleButtonDown();
    }

    void OnRightButtonDown() override
    {
      // Forward events
      vtkInteractorStyleMultiTouchCamera::OnRightButtonDown();
    }
};

#endif  // G4VTKINTERACTORSTYLE_HH
