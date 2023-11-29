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

#ifndef G4VTKGRIDPIPELINE_HH
#define G4VTKGRIDPIPELINE_HH

class G4VtkStructuredGridPipeline
{
    G4VtkStructuredGridPipeline(G4int nxIn, G4int nyIn, G4int nzIn) : nx(nxIn), ny(nyIn), nz(nzIn)
    {
      structuredGrid->SetDimensions(nx, ny, nz);
      structuredGrid->SetPoints(points);
      structuredGrid->GetCellData()->SetScalars(cellValues);
      structuredGrid->GetPointData()->SetScalars(pointValues);

      mapper->SetInputData(structuredGrid);
      actor->SetMapper(mapper);
    }

    void Modified(){};
    void Clear(){};
    void Print(){};

    G4int nx, ny, nz;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkDoubleArray> pointValues;
    vtkSmartPointer<vtkDoubleArray> cellValues;
    vtkSmartPointer<vtkStructuredGrid> structuredGrid;
    vtkSmartPointer<vtkDataSetMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
};

#endif  // G4VTKGRIDPIPELINE_HH
