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
#include "DicomBeam.hh"
#include "DicomVBeamDevice.hh"
#include "DicomBeamControlPoint.hh"
#include "DicomBeamCompensator.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeam::DicomBeam()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeam::Print( std::ostream& )
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeam::SetControlPointMetersets()
{
  G4double prevCumulativeMS = 0.;
  for( size_t ii = 0; ii < theControlPoints.size(); ii++ ){
    G4double cumulativeMS = theControlPoints[ii]->GetCumulativeMetersetWeight();
    theControlPoints[ii]->SetMetersetWeight( (cumulativeMS - prevCumulativeMS )*theMeterset);
    prevCumulativeMS = cumulativeMS;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeam::DumpToFile()
{
  std::ofstream fout("RTPlan_"+std::to_string(theNumber));
  fout << ":P BeamMeterset " << theMeterset << G4endl;

  for( size_t ii = 0; ii < theDevices.size(); ii++) {
    theDevices[ii]->DumpToFile( fout );
  }

  for( size_t ii = 0; ii < theCompensators.size(); ii++) {
    theCompensators[ii]->DumpToFile( fout );
  }

  for( size_t kk = 0; kk < theControlPoints.size(); kk++ ){
    std::string kkstr = std::to_string(theControlPoints[kk]->GetIndex());
    std::ofstream fout2("RTPlanControlPoint_"+std::to_string(theNumber)+"_"+kkstr);
    theControlPoints[kk]->DumpToFile( fout2 );
  }

}


