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
// Created by Stewart Boogert on 05/03/2023.
//

#ifndef G4VTKTEXTPIPELINE_HH
#define G4VTKTEXTPIPELINE_HH

#include "G4VVtkPipeline.hh"

class G4Text;
class G4VtkVisContext;
class vtkBillboardTextActor3D;

class G4VtkTextPipeline : public G4VVtkPipeline
{
  public:
    G4VtkTextPipeline(const G4Text& text, const G4VtkVisContext& vc,
                      const G4VisAttributes* pVisAttributes);
    ~G4VtkTextPipeline() override = default;

    void Enable() override;
    void Disable() override;

    void Print() override;
    void Modified() override;
    void Clear() override;

    virtual vtkSmartPointer<vtkBillboardTextActor3D> GetActor() { return actor; }

    virtual void SetText(const G4String& text);

    static std::size_t MakeHash(const G4Text& text, const G4VtkVisContext& vc,
                                const G4VisAttributes* pVA);

  protected:
    vtkSmartPointer<vtkBillboardTextActor3D> actor;
};

#endif  // G4VTKTEXTPIPELINE_HH
