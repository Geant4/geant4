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

#ifndef G4ChannelingFastSimInterpolation_h
#define G4ChannelingFastSimInterpolation_h

#include "G4PhysicsVector.hh"
#include "G4Physics2DVector.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicsLinearVector.hh"

/** \file G4ChannelingFastSimInterpolation.hh
* \brief Definition of the G4ChannelingFastSimInterpolation class
* The class includes spline interpolation coefficients for the important functions
* needed in Channeling FastSimulation model, i.e. electric fields, nuclear and
* electron densities and minimum energy of ionization. All the functions have
* transverse coordinates as arguments (x for crystal planes, x and y for crystal axes).
* All the functions are calculated in the co-rotating reference system (along crystal
* planes/axes).
*/

class G4ChannelingFastSimInterpolation {

public:
    G4ChannelingFastSimInterpolation(G4double dx0,
                                   G4double dy0,
                                   G4int nPointsx0,
                                   G4int nPointsy0,
                                   G4int iModel0);
    ~G4ChannelingFastSimInterpolation() = default;

    ///Get Spline Function
    G4double GetIF(G4double xx, G4double yy);

    ///Set spline coefficients
    void SetCoefficients1D(G4double AI0,
                           G4double BI0,
                           G4double CI0,
                           G4double DI0,
                           G4int i);
    void SetCoefficients2D(G4double AI3D0,
                           G4double BI3D0,
                           G4double CI3D0,
                           G4int i,
                           G4int j,
                           G4int k);

private:
    G4double Spline1D(G4double xx);
    G4double Spline2D(G4double xx, G4double yy);// cubic spline of 2-variable function

    G4double fDx=0., fDy=0.; //channel width and height
    G4double fStepi=0., fStepj=0.; //interpolation steps in x and y, respectively
    G4double fStepi2=0.; //=fStepi*fStepi
    G4int nPointsx=0, nPointsy=0; //number of interpolation nodes in x and y, respectively

    std::vector <G4double> fAI;
    std::vector <G4double> fBI;
    std::vector <G4double> fCI;
    std::vector <G4double> fDI;

    std::vector<std::vector<G4double>> fAI3D;
    std::vector<std::vector<G4double>> fBI3D;
    std::vector<std::vector<G4double>> fCI3D;
    std::vector<std::vector<G4double>> fAI3D3;
    std::vector<std::vector<G4double>> fBI3D3;
    std::vector<std::vector<G4double>> fCI3D3;

    G4int iModel=1;

};

#endif
