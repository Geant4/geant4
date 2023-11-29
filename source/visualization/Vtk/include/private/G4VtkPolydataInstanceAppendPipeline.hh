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

#ifndef G4VTKPOLYDATAINSTANCEAPPENDPIPELINE_HH
#define G4VTKPOLYDATAINSTANCEAPPENDPIPELINE_HH

#include "G4VtkPolydataInstancePipeline.hh"

class G4String;
class G4VtkVisContext;

class vtkTransformPolyDataFilter;
class vtkAppendPolyData;

class G4VtkPolydataInstanceAppendPipeline : public G4VtkPolydataInstancePipeline
{
  public:
    G4VtkPolydataInstanceAppendPipeline(G4String name, const G4VtkVisContext& vc);
    ~G4VtkPolydataInstanceAppendPipeline() override = default;

    void Print() override;

    void addInstance(G4double dx, G4double dy, G4double dz, G4double r00, G4double r01,
                     G4double r02, G4double r10, G4double r11, G4double r12, G4double r20,
                     G4double r21, G4double r22, G4double r, G4double g, G4double b, G4double a,
                     const G4String& name) override;
    void removeInstance(const G4String& name) override;

    static std::size_t MakeHash(const G4Polyhedron& p, const G4VtkVisContext& vc);

  protected:
    std::map<G4String, vtkSmartPointer<vtkTransformPolyDataFilter>> transformFilterMap;
    vtkSmartPointer<vtkAppendPolyData> appendFilter;
};

#endif  // G4VTKPOLYDATAINSTANCEAPPENDPIPELINE_HH
