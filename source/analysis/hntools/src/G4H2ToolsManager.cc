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

// Author: Ivana Hrivnacova, 18/06/2013 (ivana@ipno.in2p3.fr)

#include "G4H2ToolsManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4BaseHistoUtilities.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h2d"

#include <fstream>

using namespace G4Analysis;
using std::to_string;

//_____________________________________________________________________________
G4H2ToolsManager::G4H2ToolsManager(const G4AnalysisManagerState& state)
 : G4VH2Manager(),
   G4THnManager<tools::histo::h2d>(state, "H2")
{}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void UpdateH2Information(G4HnInformation* hnInformation,
                          const G4String& xunitName,
                          const G4String& yunitName,
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme)
{
  hnInformation->SetDimension(kX, xunitName, xfcnName, xbinScheme);
  hnInformation->SetDimension(kY, yunitName, yfcnName, ybinScheme);
}

//_____________________________________________________________________________
void AddH2Annotation(tools::histo::h2d* h2d,
                     const G4String& xunitName,
                     const G4String& yunitName,
                     const G4String& xfcnName,
                     const G4String& yfcnName)
{
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);
  UpdateTitle(yaxisTitle, yunitName, yfcnName);
  h2d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
}

//_____________________________________________________________________________
tools::histo::h2d* CreateToolsH2(
                         const G4String& title,
                         G4int nxbins, G4double xmin, G4double xmax,
                         G4int nybins, G4double ymin, G4double ymax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& xfcnName,
                         const G4String& yfcnName,
                         const G4String& xbinSchemeName,
                         const G4String& ybinSchemeName,
                         std::string_view className)

{
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);

  if ( xbinScheme != G4BinScheme::kLog && ybinScheme !=  G4BinScheme::kLog) {
    if ( xbinScheme == G4BinScheme::kUser || ybinScheme == G4BinScheme::kUser) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      Warn("User binning scheme setting was ignored.\n"
           "Linear binning will be applied with given (nbins, xmin, xmax) values.",
           className, "CreateToolsH2");
    }
    return new tools::histo::h2d(title,
                                 nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit),
                                 nybins, yfcn(ymin/yunit), yfcn(ymax/yunit));
               // h2 objects are deleted in destructor and reset when
               // closing a file.
  }
  else {
    // Compute edges
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    return new tools::histo::h2d(title, xedges, yedges);
  }
}

//_____________________________________________________________________________
tools::histo::h2d* CreateToolsH2(
                         const G4String& title,
                         const std::vector<G4double>& xedges,
                         const std::vector<G4double>& yedges,
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
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);

  return new tools::histo::h2d(title, xnewEdges, ynewEdges);
             // h2 objects are deleted in destructor and reset when
             // closing a file.
}

//_____________________________________________________________________________
void  ConfigureToolsH2(tools::histo::h2d* h2d,
                       G4int nxbins, G4double xmin, G4double xmax,
                       G4int nybins, G4double ymin, G4double ymax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName,
                       const G4String& xbinSchemeName,
                       const G4String& ybinSchemeName,
                       std::string_view className)
{
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);

  if ( xbinScheme != G4BinScheme::kLog && ybinScheme !=  G4BinScheme::kLog) {
    if ( xbinScheme == G4BinScheme::kUser || ybinScheme == G4BinScheme::kUser) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      Warn("User binning scheme setting was ignored.\n"
           "Linear binning will be applied with given (nbins, xmin, xmax) values.",
           className, "ConfigureToolsH2");
    }
    h2d->configure(nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit),
                   nybins, yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
  else {
    // Compute bins
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    h2d->configure(xedges, yedges);
  }
}

//_____________________________________________________________________________
void  ConfigureToolsH2(tools::histo::h2d* h2d,
                       const std::vector<G4double>& xedges,
                       const std::vector<G4double>& yedges,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName)
{
  auto xunit = GetUnitValue(xunitName);
  auto xfcn = GetFunction(xfcnName);
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);

  auto yunit = GetUnitValue(yunitName);
  auto yfcn = GetFunction(yfcnName);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);

  h2d->configure(xnewEdges, ynewEdges);
}

}


//
// private methods
//

//_____________________________________________________________________________
void G4H2ToolsManager::AddH2Information(const G4String& name,
                          const G4String& xunitName,
                          const G4String& yunitName,
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme) const
{
  auto hnInformation = fHnManager->AddHnInformation(name, fkDimension);
  hnInformation->AddDimension(xunitName, xfcnName, xbinScheme);
  hnInformation->AddDimension(yunitName, yfcnName, ybinScheme);
}

//
// protected methods
//

//_____________________________________________________________________________
G4int G4H2ToolsManager::CreateH2(const G4String& name,  const G4String& title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& xbinSchemeName,
                          const G4String& ybinSchemeName)

{
  Message(kVL4, "create", "H2", name);

  tools::histo::h2d* h2d
    = CreateToolsH2(title, nxbins, xmin, xmax, nybins, ymin, ymax,
                    xunitName, yunitName, xfcnName, yfcnName,
                    xbinSchemeName, ybinSchemeName, fkClass);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);

  // Save H2 information
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  AddH2Information(
    name, xunitName, yunitName, xfcnName, yfcnName, xbinScheme, ybinScheme);

  // Register histogram
  G4int id = RegisterT(h2d, name);

  Message(kVL2, "create", "H2", name);

  return id;
}

//_____________________________________________________________________________
G4int G4H2ToolsManager::CreateH2(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& xedges,
                          const std::vector<G4double>& yedges,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName)

{
  Message(kVL4, "create", "H2", name);

  tools::histo::h2d* h2d
    = CreateToolsH2(title, xedges, yedges,
                    xunitName, yunitName, xfcnName, yfcnName);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);

  // Save H2 information
  AddH2Information(
    name, xunitName, yunitName, xfcnName, yfcnName, G4BinScheme::kUser, G4BinScheme::kUser);

  // Register histogram
  G4int id = RegisterT(h2d, name);

  Message(kVL2, "create", "H2", name);

  return id;
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2(G4int id,
                            G4int nxbins, G4double xmin, G4double xmax,
                            G4int nybins, G4double ymin, G4double ymax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& xbinSchemeName,
                            const G4String& ybinSchemeName)
{
  auto h2d = GetTInFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  auto info = fHnManager->GetHnInformation(id, "SetH2");

  Message(kVL4, "configure", "H2", info->GetName());

  // Configure tools h2
  ConfigureToolsH2(
    h2d, nxbins, xmin, xmax, nybins, ymin, ymax,
    xunitName, yunitName, xfcnName, yfcnName, xbinSchemeName, ybinSchemeName,
    fkClass);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);

  // Update information
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  UpdateH2Information(
    info, xunitName, yunitName, xfcnName, yfcnName, xbinScheme, ybinScheme);

  // Set activation
  fHnManager->SetActivation(id, true);

  return true;
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2(G4int id,
                            const std::vector<G4double>& xedges,
                            const std::vector<G4double>& yedges,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName)
{
  auto h2d = GetTInFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  auto info = fHnManager->GetHnInformation(id, "SetH2");

  Message(kVL4, "configure", "H2", info->GetName());

  // Configure tools h2
  ConfigureToolsH2(h2d, xedges, yedges, xunitName, yunitName, xfcnName, yfcnName);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);

  // Update information
  UpdateH2Information(
    info, xunitName, yunitName, xfcnName, yfcnName, G4BinScheme::kUser, G4BinScheme::kUser);

  // Set activation
  fHnManager->SetActivation(id, true);

  return true;
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::ScaleH2(G4int id, G4double factor)
{
  auto h2d = GetTInFunction(id, "ScaleH2", false, false);
  if ( ! h2d ) return false;

  return h2d->scale(factor);
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::FillH2(G4int id,
                                G4double xvalue, G4double yvalue,
                                G4double weight)
{
  auto h2d = GetTInFunction(id, "FillH2", true, false);
  if ( ! h2d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    return false;
  }

  G4HnDimensionInformation* xInfo
    = fHnManager->GetHnDimensionInformation(id, kX, "FillH2");
  G4HnDimensionInformation* yInfo
    = fHnManager->GetHnDimensionInformation(id, kY, "FillH2");

  h2d->fill(xInfo->fFcn(xvalue/xInfo->fUnit),
            yInfo->fFcn(yvalue/yInfo->fUnit), weight);

  if ( IsVerbose(kVL4) ) {
    Message(kVL4, "fill", "H2",
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
G4int  G4H2ToolsManager::GetH2Id(const G4String& name, G4bool warn) const
{
  return GetTId(name, warn);
}

//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nxbins(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2NXbins");
  if ( ! h2d ) return 0;

  return GetNbins(*h2d, kX);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h2d = GetTInFunction(id, "GetH2Xmin");
  if ( ! h2d ) return 0.;

  return GetMin(*h2d, kX);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmax(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2Xmax");
  if ( ! h2d ) return 0.;

  return GetMax(*h2d, kX);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2XWidth(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2XWidth", true, false);
  if ( ! h2d ) return 0.;

  return GetWidth(*h2d, kX, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nybins(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2NYbins");
  if ( ! h2d ) return 0;

  return GetNbins(*h2d, kY);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h2d = GetTInFunction(id, "GetH2Ymin");
  if ( ! h2d ) return 0.;

  return GetMin(*h2d, kY);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymax(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2Ymax");
  if ( ! h2d ) return 0.;

  return GetMax(*h2d, kY);
}

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2YWidth(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2YWidth", true, false);
  if ( ! h2d ) return 0.;

  return GetWidth(*h2d, kY, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2Title(G4int id, const G4String& title)
{
  auto h2d = GetTInFunction(id, "SetH2Title");
  if ( ! h2d ) return false;

  return SetTitle(*h2d, title);
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  auto h2d = GetTInFunction(id, "SetH2XAxisTitle");
  if ( ! h2d ) return false;

  return SetAxisTitle(*h2d, kX, title);
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  auto h2d = GetTInFunction(id, "SetH2YAxisTitle");
  if ( ! h2d ) return false;

  return SetAxisTitle(*h2d, kY, title);
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  auto h2d = GetTInFunction(id, "SetH2ZAxisTitle");
  if ( ! h2d ) return false;

  return SetAxisTitle(*h2d, kZ, title);
}

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2Title(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2Title");
  if ( ! h2d ) return "";

  return GetTitle(*h2d);
}

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2XAxisTitle(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2XAxisTitle");
  if ( ! h2d ) return "";

  return GetAxisTitle(*h2d, kX, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2YAxisTitle(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2YAxisTitle");
  if ( ! h2d ) return "";

  return GetAxisTitle(*h2d, kY, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2ZAxisTitle(G4int id) const
{
  auto h2d = GetTInFunction(id, "GetH2ZAxisTitle");
  if ( ! h2d ) return "";

  return GetAxisTitle(*h2d, kZ, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4bool G4H2ToolsManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file

  // Do nothing if no histograms are selected
  if ( ! fHnManager->IsAscii() ) return true;

  // Write h2 histograms
  for ( G4int i=0; i<G4int(fTVector.size()); ++i ) {
    auto id = i + fHnManager->GetFirstId();
    auto info = fHnManager->GetHnInformation(id,"WriteOnAscii");
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->GetAscii() ) continue;
    auto h2 = fTVector[i];

    Message(kVL3, "write on ascii", "h2d", info->GetName());

    output << "\n  2D histogram " << id << ": " << h2->title()
           << "\n \n \t \t     X \t\t     Y \t\t Bin Height" << G4endl;

    for (G4int j=0; j< G4int(h2->axis_x().bins()); ++j) {
      for (G4int k=0; k< G4int(h2->axis_y().bins()); ++k) {
        output << "  " << j << "\t" << k << "\t"
               << h2->axis_x().bin_center(j) << "\t"
               << h2->axis_y().bin_center(k) << "\t"
               << h2->bin_height(j, k) << G4endl;
      }
    }
  }

  return output.good();
}

//
// public methods
//

//_____________________________________________________________________________
G4int G4H2ToolsManager::AddH2(const G4String& name, tools::histo::h2d* h2d)
{
  Message(kVL4, "add", "H2", name);

  // Add annotation
  AddH2Annotation(h2d, "none", "none", "none", "none");
  // Add information
  AddH2Information(name, "none", "none", "none", "none",
                   G4BinScheme::kLinear, G4BinScheme::kLinear);

  // Register histogram
  G4int id = RegisterT(h2d, name);

  Message(kVL2, "add", "H2", name);

  return id;
}

//_____________________________________________________________________________
void G4H2ToolsManager::AddH2Vector(
                          const std::vector<tools::histo::h2d*>& h2Vector)
{
  AddTVector(h2Vector);
}
//_____________________________________________________________________________
tools::histo::h2d*  G4H2ToolsManager::GetH2(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const
{
  return GetTInFunction(id, "GetH2", warn, onlyIfActive);
}

