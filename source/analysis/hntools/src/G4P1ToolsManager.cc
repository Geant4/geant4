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

// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#include "G4P1ToolsManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4BaseHistoUtilities.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/p1d"

#include <fstream>

using namespace G4Analysis;
using std::to_string;

//
// Constructors, destructor
//

//_____________________________________________________________________________
G4P1ToolsManager::G4P1ToolsManager(const G4AnalysisManagerState& state)
 : G4VP1Manager(),
   G4THnManager<tools::histo::p1d>(state, "P1")
{}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void UpdateP1Information(G4HnInformation* hnInformation,
                          const G4String& xunitName,
                          const G4String& yunitName,
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme)
{
  hnInformation->SetDimension(kX, xunitName, xfcnName, xbinScheme);
  hnInformation->SetDimension(kY, yunitName, yfcnName, G4BinScheme::kLinear);
}

//_____________________________________________________________________________
void AddP1Annotation(tools::histo::p1d* p1d,
                     const G4String& xunitName,
                     const G4String& yunitName,
                     const G4String& xfcnName,
                     const G4String& yfcnName)
{
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);
  UpdateTitle(yaxisTitle, yunitName, yfcnName);
  p1d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  p1d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
}

//_____________________________________________________________________________
tools::histo::p1d* CreateToolsP1(const G4String& title,
                         G4int nbins, G4double xmin, G4double xmax,
                         G4double ymin, G4double ymax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& xfcnName,
                         const G4String& yfcnName,
                         const G4String& xbinSchemeName,
                         std::string_view className)
{
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);

  if ( xbinScheme != G4BinScheme::kLog ) {
    if ( xbinScheme == G4BinScheme::kUser ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      Warn("User binning scheme setting was ignored.\n"
           "Linear binning will be applied with given (nbins, xmin, xmax) values.",
           className, "CreateToolsP1");
    }
    if ( ymin == 0. && ymax == 0.) {
      return new tools::histo::p1d(title,
                                 nbins, xfcn(xmin/xunit), xfcn(xmax/xunit));
    } else {
      return new tools::histo::p1d(title,
                                 nbins, xfcn(xmin/xunit), xfcn(xmax/xunit),
                                 yfcn(ymin/yunit), yfcn(ymax/yunit));
    }
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    if ( ymin == 0. && ymax == 0.) {
      return new tools::histo::p1d(title, edges);
    } else {
      return new tools::histo::p1d(title, edges, yfcn(ymin/yunit), yfcn(ymax/yunit));
    }
  }
}

//_____________________________________________________________________________
tools::histo::p1d* CreateToolsP1(const G4String& title,
                         const std::vector<G4double>& edges,
                         G4double ymin, G4double ymax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& xfcnName,
                         const G4String& yfcnName)
{
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);

  // Apply function
  std::vector<G4double> newEdges;
  ComputeEdges(edges, xunit, xfcn, newEdges);

  return new tools::histo::p1d(title, newEdges, yfcn(ymin/yunit), yfcn(ymax/yunit));
}

//_____________________________________________________________________________
void ConfigureToolsP1(tools::histo::p1d* p1d,
                       G4int nbins, G4double xmin, G4double xmax,
                       G4double ymin, G4double ymax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName,
                       const G4String& xbinSchemeName,
                       std::string_view className)
{
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);

  if ( xbinScheme != G4BinScheme::kLog ) {
    if ( xbinScheme == G4BinScheme::kUser ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      Warn("User binning scheme setting was ignored.\n"
           "Linear binning will be applied with given (nbins, xmin, xmax) values.",
           className, "ConfigureToolsP1");
    }
    if ( ymin == 0. && ymax == 0. ) {
      p1d->configure(nbins, xfcn(xmin/xunit), xfcn(xmax/xunit));
    } else {
      p1d->configure(nbins, xfcn(xmin/xunit), xfcn(xmax/xunit),
                     yfcn(ymin/yunit), yfcn(ymax/yunit));
    }
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    if ( ymin == 0. && ymax == 0. ) {
      p1d->configure(edges);
    } else {
      p1d->configure(edges, yfcn(ymin/yunit), yfcn(ymax/yunit));
    }
  }
}

//_____________________________________________________________________________
void ConfigureToolsP1(tools::histo::p1d* p1d,
                      const std::vector<G4double>& edges,
                      G4double ymin, G4double ymax,
                      const G4String& xunitName,
                      const G4String& yunitName,
                      const G4String& xfcnName,
                      const G4String& yfcnName)
{
  // Apply function to edges
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, xunit, xfcn, newEdges);

  if ( ymin == 0. && ymax == 0. ) {
    p1d->configure(newEdges);
  } else {
    p1d->configure(newEdges, yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
}

}

//
// private methods
//

//_____________________________________________________________________________
void G4P1ToolsManager::AddP1Information(const G4String& name,
                          const G4String& xunitName,
                          const G4String& yunitName,
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme) const
{
  auto hnInformation = fHnManager->AddHnInformation(name, fkDimension);
  hnInformation->AddDimension(xunitName, xfcnName, xbinScheme);
  hnInformation->AddDimension(yunitName, yfcnName, G4BinScheme::kLinear);
}

//
// protected methods
//

//_____________________________________________________________________________
G4int G4P1ToolsManager::CreateP1(const G4String& name,  const G4String& title,
                          G4int nbins, G4double xmin, G4double xmax,
                          G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& xbinSchemeName)
{
  Message(kVL4, "create", "P1", name);

  tools::histo::p1d* p1d
    = CreateToolsP1(title, nbins, xmin, xmax, ymin, ymax,
                    xunitName, yunitName, xfcnName, yfcnName,
                    xbinSchemeName, fkClass);

  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);

  // Save P1 information
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  AddP1Information(
    name, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Register profile
  G4int id = RegisterT(p1d, name);

  Message(kVL2, "create", "P1", name);

  return id;
}

//_____________________________________________________________________________
G4int G4P1ToolsManager::CreateP1(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& edges,
                          G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName)
{
  Message(kVL4, "create", "P1", name);

  tools::histo::p1d* p1d
    = CreateToolsP1(title, edges, ymin, ymax,
                    xunitName, yunitName, xfcnName, yfcnName);

  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);

  // Save P1 information
  AddP1Information(
    name, xunitName, yunitName, xfcnName, yfcnName, G4BinScheme::kUser);

  // Register profile
  G4int id = RegisterT(p1d, name);

  Message(kVL2, "create", "P1", name);

  return id;
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1(G4int id,
                            G4int nbins, G4double xmin, G4double xmax,
                            G4double ymin, G4double ymax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& xbinSchemeName)
{
  auto p1d = GetTInFunction(id, "SetP1", false, false);
  if ( ! p1d ) return false;

  auto info = fHnManager->GetHnInformation(id,"SetP1");

  Message(kVL4, "configure", "P1", info->GetName());

  // Configure tools p1
  ConfigureToolsP1(
    p1d, nbins, xmin, xmax, ymin, ymax,
    xunitName, yunitName, xfcnName, yfcnName, xbinSchemeName, fkClass);

  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);

  // Update information
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  UpdateP1Information(
    info, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Set activation
  fHnManager->SetActivation(id, true);

  return true;
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1(G4int id,
                           const std::vector<G4double>& edges,
                           G4double ymin, G4double ymax,
                           const G4String& xunitName, const G4String& yunitName,
                           const G4String& xfcnName, const G4String& yfcnName)
{
  auto p1d = GetTInFunction(id, "SetP1", false, false);
  if ( ! p1d ) return false;

  auto info = fHnManager->GetHnInformation(id,"SetP1");

  Message(kVL4, "configure", "P1", info->GetName());

  // Configure tools p1
  ConfigureToolsP1(p1d, edges, ymin, ymax,
                   xunitName, yunitName, xfcnName, yfcnName);

  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);

  // Update information
  UpdateP1Information(
    info, xunitName, yunitName, xfcnName, yfcnName, G4BinScheme::kUser);

  // Set activation
  fHnManager->SetActivation(id, true);

  return true;
}


//_____________________________________________________________________________
G4bool G4P1ToolsManager::ScaleP1(G4int id, G4double factor)
{
  auto p1d = GetTInFunction(id, "ScaleP1", false, false);
  if ( ! p1d ) return false;

  return p1d->scale(factor);
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::FillP1(G4int id, G4double xvalue, G4double yvalue,
                                G4double weight)
{
  auto p1d = GetTInFunction(id, "FillP1", true, false);
  if ( ! p1d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillP1 for " << id << G4endl;
    return false;
  }

  auto xInfo
    = fHnManager->GetHnDimensionInformation(id, kX, "FillP1");
  auto yInfo
    = fHnManager->GetHnDimensionInformation(id, kY, "FillP1");

  p1d->fill(xInfo->fFcn(xvalue/xInfo->fUnit),
            yInfo->fFcn(yvalue/yInfo->fUnit), weight);

  if ( IsVerbose(kVL4) ) {
    Message(kVL4, "fill", "P1",
      " id " + to_string(id) +
      " xvalue " + to_string(xvalue) +
      " xfcn(xvalue/xunit) " +  to_string(xInfo->fFcn(xvalue/xInfo->fUnit)) +
      " yvalue " + to_string(yvalue) +
      " yfcn(yvalue/yunit) " +  to_string(yInfo->fFcn(yvalue/yInfo->fUnit)) +
      " weight " + to_string(weight));
  }

  return true;
}

//_____________________________________________________________________________
G4int  G4P1ToolsManager::GetP1Id(const G4String& name, G4bool warn) const
{
  return GetTId(name, warn);
}

//_____________________________________________________________________________
G4int G4P1ToolsManager::GetP1Nbins(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1Nbins");
  if ( ! p1d ) return 0;

  return GetNbins(*p1d, kX);
}

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Xmin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  auto p1d = GetTInFunction(id, "GetP1Xmin");
  if ( ! p1d ) return 0.;

  return GetMin(*p1d, kX);
}

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Xmax(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1Xmax");
  if ( ! p1d ) return 0.;

  return GetMax(*p1d, kX);
}

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1XWidth(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1XWidth", true, false);
  if ( ! p1d ) return 0.;

  return GetWidth(*p1d, kX, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Ymin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  auto p1d = GetTInFunction(id, "GetP1Ymin");
  if ( ! p1d ) return 0.;

  return p1d->min_v();
}

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Ymax(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1Ymax");
  if ( ! p1d ) return 0.;

  return p1d->max_v();
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1Title(G4int id, const G4String& title)
{
  auto p1d = GetTInFunction(id, "SetP1Title");
  if ( ! p1d ) return false;

  return SetTitle(*p1d, title);
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1XAxisTitle(G4int id, const G4String& title)
{
  auto p1d = GetTInFunction(id, "SetP1XAxisTitle");
  if ( ! p1d ) return false;

  return SetAxisTitle(*p1d, kX, title);
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1YAxisTitle(G4int id, const G4String& title)
{
  auto p1d = GetTInFunction(id, "SetP1YAxisTitle");
  if ( ! p1d ) return false;

  return SetAxisTitle(*p1d, kY, title);
}

//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1Title(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1Title");
  if ( ! p1d ) return "";

  return GetTitle(*p1d);
}


//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1XAxisTitle(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1XAxisTitle");
  if ( ! p1d ) return "";

  return GetAxisTitle(*p1d, kX, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1YAxisTitle(G4int id) const
{
  auto p1d = GetTInFunction(id, "GetP1YAxisTitle");
  if ( ! p1d ) return "";

  return GetAxisTitle(*p1d, kY, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file

  // Do nothing if no histograms are selected
  if ( ! fHnManager->IsAscii() ) return true;

  // Write p1 histograms
  for ( G4int i=0; i<G4int(fTVector.size()); ++i ) {
    auto id = i + fHnManager->GetFirstId();
    auto info = fHnManager->GetHnInformation(id,"WriteOnAscii");
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->GetAscii() ) continue;
    auto p1 = fTVector[i];

    Message(kVL3, "write on ascii", "p1d", info->GetName());

    output << "\n  1D profile " << id << ": " << p1->title()
           << "\n \n \t     X \t\t    MeanY" << G4endl;

    for (G4int j=0; j< G4int(p1->axis().bins()); ++j) {
       auto sw = p1->bin_Sw(j);
       auto svw = p1->bin_Svw(j);
       auto mean = ( sw != 0. ) ?  (svw / sw) : 0.;
       output << "  " << j << "\t"
              << p1->axis().bin_center(j) << "\t"
              << mean << G4endl;
    }
  }

  return output.good();
}

//
// public methods
//

//_____________________________________________________________________________
G4int G4P1ToolsManager::AddP1(const G4String& name, tools::histo::p1d* p1d)
{
  Message(kVL4, "add", "P1", name);

  // Add annotation
  AddP1Annotation(p1d, "none", "none", "none", "none");
  // Add information
  AddP1Information(name, "none", "none", "none", "none", G4BinScheme::kLinear);

  // Register profile
  auto id = RegisterT(p1d, name);

  Message(kVL2, "add", "P1", name);

  return id;
}

//_____________________________________________________________________________
void G4P1ToolsManager::AddP1Vector(
                          const std::vector<tools::histo::p1d*>& p1Vector)
{
  AddTVector(p1Vector);
}

//_____________________________________________________________________________
tools::histo::p1d*  G4P1ToolsManager::GetP1(G4int id, G4bool warn,
                                            G4bool onlyIfActive) const
{
  return GetTInFunction(id, "GetP1", warn, onlyIfActive);
}

