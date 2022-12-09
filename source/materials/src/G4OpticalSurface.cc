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
//
////////////////////////////////////////////////////////////////////////
// Optical Surface Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpticalSurface.cc
// Description: An optical surface class for use in G4OpBoundaryProcess
// Version:     2.0
// Created:     1997-06-26
// Author:      Peter Gumplinger
// updated:     2017-02-24 Mariele Stockhoff add DAVIS model
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <zlib.h>

#include "globals.hh"
#include "G4OpticalSurface.hh"

G4OpticalSurface& G4OpticalSurface::operator=(const G4OpticalSurface& right)
{
  if(this != &right)
  {
    theName                    = right.theName;
    theType                    = right.theType;
    theModel                   = right.theModel;
    theFinish                  = right.theFinish;
    sigma_alpha                = right.sigma_alpha;
    polish                     = right.polish;
    theMaterialPropertiesTable = right.theMaterialPropertiesTable;

    delete[] AngularDistribution;
    AngularDistribution =
      new G4float[incidentIndexMax * thetaIndexMax * phiIndexMax];
    *(AngularDistribution) = *(right.AngularDistribution);

    delete[] AngularDistributionLUT;
    AngularDistributionLUT    = new G4float[indexmax];
    *(AngularDistributionLUT) = *(right.AngularDistributionLUT);

    delete[] Reflectivity;
    Reflectivity    = new G4float[RefMax];
    *(Reflectivity) = *(right.Reflectivity);

    delete DichroicVector;
    DichroicVector  = new G4Physics2DVector();
    *DichroicVector = *(right.DichroicVector);
  }
  return *this;
}

G4OpticalSurface::G4OpticalSurface(const G4String& name,
                                   G4OpticalSurfaceModel model,
                                   G4OpticalSurfaceFinish finish,
                                   G4SurfaceType type, G4double value)
  : G4SurfaceProperty(name, type)
  , theModel(model)
  , theFinish(finish)
  , theMaterialPropertiesTable(nullptr)
{
  AngularDistribution = nullptr;

  AngularDistributionLUT = nullptr;
  Reflectivity           = nullptr;

  DichroicVector = nullptr;

  switch(theModel)
  {
    case glisur:
      polish      = value;
      sigma_alpha = 0.0;
      break;
    case LUT:
    case dichroic:
    case DAVIS:
      ReadDataFile();
      // fall through
    case unified:
      sigma_alpha = value;
      polish      = 0.0;
      break;
    default:
      G4Exception("G4OpticalSurface::G4OpticalSurface()", "mat309",
                  FatalException, "Constructor called with INVALID model.");
  }
}

G4OpticalSurface::~G4OpticalSurface()
{
    delete[] AngularDistribution;

    delete[] AngularDistributionLUT;

    delete[] Reflectivity;

    delete DichroicVector;
}

G4OpticalSurface::G4OpticalSurface(const G4OpticalSurface& right)
  : G4SurfaceProperty(right.theName, right.theType)
{
  *this                            = right;
  this->theName                    = right.theName;
  this->theType                    = right.theType;
  this->theModel                   = right.theModel;
  this->theFinish                  = right.theFinish;
  this->sigma_alpha                = right.sigma_alpha;
  this->polish                     = right.polish;
  this->theMaterialPropertiesTable = right.theMaterialPropertiesTable;

    delete[] AngularDistribution;
  this->AngularDistribution =
    new G4float[incidentIndexMax * thetaIndexMax * phiIndexMax];
  *(this->AngularDistribution) = *(right.AngularDistribution);

    delete[] AngularDistributionLUT;
  this->AngularDistributionLUT    = new G4float[indexmax];
  *(this->AngularDistributionLUT) = *(right.AngularDistributionLUT);

    delete[] Reflectivity;
  this->Reflectivity    = new G4float[RefMax];
  *(this->Reflectivity) = *(right.Reflectivity);

    delete DichroicVector;
  this->DichroicVector    = new G4Physics2DVector();
  *(this->DichroicVector) = *(right.DichroicVector);
}

G4bool G4OpticalSurface::operator==(const G4OpticalSurface& right) const
{
  return (this == (G4OpticalSurface*) &right);
}

G4bool G4OpticalSurface::operator!=(const G4OpticalSurface& right) const
{
  return (this != (G4OpticalSurface*) &right);
}

G4int G4OpticalSurface::GetInmax() const { return indexmax; }

G4int G4OpticalSurface::GetLUTbins() const { return LUTbins; }

G4int G4OpticalSurface::GetRefMax() const { return RefMax; }

G4int G4OpticalSurface::GetThetaIndexMax() const { return thetaIndexMax; }

G4int G4OpticalSurface::GetPhiIndexMax() const { return phiIndexMax; }

void G4OpticalSurface::DumpInfo() const
{
  // Dump info for surface

  G4cout << "  Surface type   = " << G4int(theType) << G4endl
         << "  Surface finish = " << G4int(theFinish) << G4endl
         << "  Surface model  = " << G4int(theModel) << G4endl << G4endl
         << "  Surface parameter " << G4endl << "  ----------------- "
         << G4endl;

  if(theModel == glisur)
  {
    G4cout << " polish: " << polish << G4endl;
  }
  else
  {
    G4cout << " sigma_alpha: " << sigma_alpha << G4endl;
  }
  G4cout << G4endl;
}

void G4OpticalSurface::SetType(const G4SurfaceType& type)
{
  theType = type;
  ReadDataFile();
}

void G4OpticalSurface::SetFinish(const G4OpticalSurfaceFinish finish)
{
  theFinish = finish;
  ReadDataFile();
}

void G4OpticalSurface::ReadDataFile()
{
  // type and finish can be set in either order. Thus, we can't check
  // for consistency. Need to read file on setting either type or finish.
  switch(theType)
  {
    case dielectric_LUT:
      if(AngularDistribution == nullptr)
      {
        AngularDistribution =
          new G4float[incidentIndexMax * thetaIndexMax * phiIndexMax];
      }
      ReadLUTFile();
      break;
    case dielectric_LUTDAVIS:
      if(AngularDistributionLUT == nullptr)
      {
        AngularDistributionLUT = new G4float[indexmax];
      }
      ReadLUTDAVISFile();

      if(Reflectivity == nullptr)
      {
        Reflectivity = new G4float[RefMax];
      }
      ReadReflectivityLUTFile();
      break;
    case dielectric_dichroic:
      if(DichroicVector == nullptr)
      {
        DichroicVector = new G4Physics2DVector();
      }
      ReadDichroicFile();
      break;
    default:
      break;
  }
}

void G4OpticalSurface::ReadLUTFile()
{
  G4String readLUTFileName;

  switch(theFinish)
  {
    case polishedlumirrorglue:
      readLUTFileName = "PolishedLumirrorGlue.z";
      break;
    case polishedlumirrorair:
      readLUTFileName = "PolishedLumirror.z";
      break;
    case polishedteflonair:
      readLUTFileName = "PolishedTeflon.z";
      break;
    case polishedtioair:
      readLUTFileName = "PolishedTiO.z";
      break;
    case polishedtyvekair:
      readLUTFileName = "PolishedTyvek.z";
      break;
    case polishedvm2000glue:
      readLUTFileName = "PolishedVM2000Glue.z";
      break;
    case polishedvm2000air:
      readLUTFileName = "PolishedVM2000.z";
      break;
    case etchedlumirrorglue:
      readLUTFileName = "EtchedLumirrorGlue.z";
      break;
    case etchedlumirrorair:
      readLUTFileName = "EtchedLumirror.z";
      break;
    case etchedteflonair:
      readLUTFileName = "EtchedTeflon.z";
      break;
    case etchedtioair:
      readLUTFileName = "EtchedTiO.z";
      break;
    case etchedtyvekair:
      readLUTFileName = "EtchedTyvek.z";
      break;
    case etchedvm2000glue:
      readLUTFileName = "EtchedVM2000Glue.z";
      break;
    case etchedvm2000air:
      readLUTFileName = "EtchedVM2000.z";
      break;
    case groundlumirrorglue:
      readLUTFileName = "GroundLumirrorGlue.z";
      break;
    case groundlumirrorair:
      readLUTFileName = "GroundLumirror.z";
      break;
    case groundteflonair:
      readLUTFileName = "GroundTeflon.z";
      break;
    case groundtioair:
      readLUTFileName = "GroundTiO.z";
      break;
    case groundtyvekair:
      readLUTFileName = "GroundTyvek.z";
      break;
    case groundvm2000glue:
      readLUTFileName = "GroundVM2000Glue.z";
      break;
    case groundvm2000air:
      readLUTFileName = "GroundVM2000.z";
      break;
    default:
      return;
  }

  std::istringstream iss;
  ReadCompressedFile(readLUTFileName, iss);

  size_t idxmax = incidentIndexMax * thetaIndexMax * phiIndexMax;
  for(size_t i = 0; i < idxmax; ++i)
  {
    iss >> AngularDistribution[i];
  }
  G4cout << "LUT - data file: " << readLUTFileName << " read in! " << G4endl;
}

void G4OpticalSurface::ReadLUTDAVISFile()
{
  G4String readLUTDAVISFileName;

  switch(theFinish)
  {
    case Rough_LUT:
      readLUTDAVISFileName = "Rough_LUT.z";
      break;
    case RoughTeflon_LUT:
      readLUTDAVISFileName = "RoughTeflon_LUT.z";
      break;
    case RoughESR_LUT:
      readLUTDAVISFileName = "RoughESR_LUT.z";
      break;
    case RoughESRGrease_LUT:
      readLUTDAVISFileName = "RoughESRGrease_LUT.z";
      break;
    case Polished_LUT:
      readLUTDAVISFileName = "Polished_LUT.z";
      break;
    case PolishedTeflon_LUT:
      readLUTDAVISFileName = "PolishedTeflon_LUT.z";
      break;
    case PolishedESR_LUT:
      readLUTDAVISFileName = "PolishedESR_LUT.z";
      break;
    case PolishedESRGrease_LUT:
      readLUTDAVISFileName = "PolishedESRGrease_LUT.z";
      break;
    case Detector_LUT:
      readLUTDAVISFileName = "Detector_LUT.z";
      break;
    default:
      return;
  }

  std::istringstream iss;
  ReadCompressedFile(readLUTDAVISFileName, iss);

  for(size_t i = 0; i < indexmax; ++i)
  {
    iss >> AngularDistributionLUT[i];
  }
  G4cout << "LUT DAVIS - data file: " << readLUTDAVISFileName << " read in! "
         << G4endl;
}

void G4OpticalSurface::ReadReflectivityLUTFile()
{
  G4String readReflectivityLUTFileName;

  switch(theFinish)
  {
    case Rough_LUT:
      readReflectivityLUTFileName = "Rough_LUTR.z";
      break;
    case RoughTeflon_LUT:
      readReflectivityLUTFileName = "RoughTeflon_LUTR.z";
      break;
    case RoughESR_LUT:
      readReflectivityLUTFileName = "RoughESR_LUTR.z";
      break;
    case RoughESRGrease_LUT:
      readReflectivityLUTFileName = "RoughESRGrease_LUTR.z";
      break;
    case Polished_LUT:
      readReflectivityLUTFileName = "Polished_LUTR.z";
      break;
    case PolishedTeflon_LUT:
      readReflectivityLUTFileName = "PolishedTeflon_LUTR.z";
      break;
    case PolishedESR_LUT:
      readReflectivityLUTFileName = "PolishedESR_LUTR.z";
      break;
    case PolishedESRGrease_LUT:
      readReflectivityLUTFileName = "PolishedESRGrease_LUTR.z";
      break;
    case Detector_LUT:
      readReflectivityLUTFileName = "Detector_LUTR.z";
      break;
    default:
      return;
  }

  std::istringstream iss;
  ReadCompressedFile(readReflectivityLUTFileName, iss);

  for(size_t i = 0; i < RefMax; ++i)
  {
    iss >> Reflectivity[i];
  }
  G4cout << "LUT DAVIS - reflectivity data file: "
         << readReflectivityLUTFileName << " read in! " << G4endl;
}

// uncompress one data file into the input string stream
void G4OpticalSurface::ReadCompressedFile(const G4String& filename,
                                          std::istringstream& iss)
{
  G4String* dataString  = nullptr;
  G4String path         = G4FindDataDir("G4REALSURFACEDATA");
  G4String compfilename = path + "/" + filename;
  // create input stream with binary mode operation and position at end of file
  std::ifstream in(compfilename, std::ios::binary | std::ios::ate);
  if(in.good())
  {
    // get current position in the stream (was set to the end)
    G4int fileSize = (G4int)in.tellg();
    // set current position being the beginning of the stream
    in.seekg(0, std::ios::beg);
    // create (zlib) byte buffer for the data
    Bytef* compdata = new Bytef[fileSize];
    while(in)
    {
      in.read((char*) compdata, fileSize);
    }
    // create (zlib) byte buffer for the uncompressed data
    uLongf complen    = (uLongf)(fileSize * 4);
    Bytef* uncompdata = new Bytef[complen];
    while(Z_OK != uncompress(uncompdata, &complen, compdata, fileSize))
    {
      // increase uncompressed byte buffer
      delete[] uncompdata;
      complen *= 2;
      uncompdata = new Bytef[complen];
    }
    // delete the compressed data buffer
    delete[] compdata;
    // create a string from uncompressed data (will be deallocated by caller)
    dataString = new G4String((char*) uncompdata, (long) complen);
    // delete the uncompressed data buffer
    delete[] uncompdata;
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "Problem while trying to read " + compfilename + " data file.\n";
    G4Exception("G4OpticalSurface::ReadCompressedFile", "mat316",
                FatalException, ed);
    return;
  }
  // create the input string stream from the data string
  if(dataString != nullptr)
  {
    iss.str(*dataString);
    in.close();
    delete dataString;
    G4cout << "G4OpticalSurface: data file " << compfilename
           << " successfully read in." << G4endl;
  }
}

void G4OpticalSurface::ReadDichroicFile()
{
  const char* datadir = G4FindDataDir("G4DICHROICDATA");

  if(datadir == nullptr)
  {
    G4Exception("G4OpticalSurface::ReadDichroicFile()", "mat313",
                FatalException,
                "Environment variable G4DICHROICDATA not defined");
    return;
  }

  std::ostringstream ost;
  ost << datadir;
  std::ifstream fin(ost.str().c_str());
  if(!fin.is_open())
  {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4OpticalSurface::ReadDichroicFile()", "mat314",
                FatalException, ed, " ");
    return;
  }

  if(!(DichroicVector->Retrieve(fin)))
  {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4OpticalSurface::ReadDichroicFile()", "mat315",
                FatalException, ed, " ");
    return;
  }

  //  DichroicVector->SetBicubicInterpolation(true);

  G4cout << " *** Dichroic surface data file *** " << G4endl;

  G4int numberOfXNodes = (G4int)DichroicVector->GetLengthX();
  G4int numberOfYNodes = (G4int)DichroicVector->GetLengthY();

  G4cout << "numberOfXNodes: " << numberOfXNodes << G4endl;
  G4cout << "numberOfYNodes: " << numberOfYNodes << G4endl;

  if(0 > numberOfXNodes || numberOfXNodes >= INT_MAX)
  {
    numberOfXNodes = 0;
  }
  if(0 > numberOfYNodes || numberOfYNodes >= INT_MAX)
  {
    numberOfYNodes = 0;
  }

  G4PV2DDataVector xVector;
  G4PV2DDataVector yVector;

  xVector.resize(numberOfXNodes, 0.);
  yVector.resize(numberOfYNodes, 0.);

  for(G4int i = 0; i < numberOfXNodes; ++i)
  {
    G4cout << "i: " << DichroicVector->GetX(i) << G4endl;
    xVector[i] = DichroicVector->GetX(i);
  }
  for(G4int j = 0; j < numberOfYNodes; ++j)
  {
    G4cout << "j: " << DichroicVector->GetY(j) << G4endl;
    yVector[j] = DichroicVector->GetY(j);
  }

  for(G4int j = 0; j < numberOfYNodes; ++j)
  {
    for(G4int i = 0; i < numberOfXNodes; ++i)
    {
      G4cout << " i: " << i << " j: " << j << " "
             << DichroicVector->GetValue(i, j) << G4endl;
    }
  }
}
