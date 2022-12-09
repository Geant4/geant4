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
// 
// John Allison  5th September 2018, based on G4PhysicalVolumeSearchScene
// An artificial scene to find physical volumes. Instead of returning the
// first occurrence (G4PhysicalVolumeSearchScene) this class (note the extra
// 's' in the name of this class) returns a vector of all occurrences.

#include "G4PhysicalVolumesSearchScene.hh"

#include <regex>

G4PhysicalVolumesSearchScene::G4PhysicalVolumesSearchScene
(G4PhysicalVolumeModel* pSearchVolumesModel,
 const G4String&        requiredPhysicalVolumeName,
 G4int                  requiredCopyNo)
: fpSearchVolumesModel  (pSearchVolumesModel)
, fMatcher              (requiredPhysicalVolumeName)
, fRequiredCopyNo       (requiredCopyNo)
{}

void G4PhysicalVolumesSearchScene::ProcessVolume (const G4VSolid&)
{
  G4VPhysicalVolume* pCurrentPV = fpSearchVolumesModel->GetCurrentPV();
  const G4String& name = pCurrentPV->GetName();
  G4int copyNo = fpSearchVolumesModel->GetCurrentPVCopyNo();

  // Match the name with the required physical volume name. The latter can be of
  // the form "/regexp/", where regexp is a regular expression (see C++ regex),
  // or a plain name, in which case there must be an exact match.
  if (fMatcher.Match(name)) {
    if ((fRequiredCopyNo < 0 ||  // I.e., ignore negative request
         fRequiredCopyNo == copyNo)) {
      auto basePath = fpSearchVolumesModel->GetFullPVPath();
      basePath.pop_back();  // Base node is one up from found node
      // Mark base path nodes as not drawn
      for (auto& node: basePath) node.SetDrawn(false);
      fFindings.push_back
      (Findings
       (fpSearchVolumesModel->GetTopPhysicalVolume(),
        pCurrentPV,
        copyNo,
        fpSearchVolumesModel->GetCurrentDepth(),
        basePath,
        fpSearchVolumesModel->GetFullPVPath(),
        *fpCurrentObjectTransformation));
    }
  }
}

G4PhysicalVolumesSearchScene::Matcher::Matcher(const G4String& requiredMatch)
: fRegexFlag(false)
{
  if (requiredMatch.size()) {
    std::size_t last = requiredMatch.size() - 1;
    // If required name begins and ends with '/', treat as a regular expression.
    // 0 causes a conversion ambiguity that upsets the Windows compiler, so use 0U.
    if (requiredMatch[0U] == '/' && requiredMatch[(G4int)last] == '/') {
      if (last > 1) {  // Non-null regexp
        // regex match required
        fRegexFlag = true;
        // Extract the required regex
        fRequiredMatch = requiredMatch.substr(1,last-1);
      }
    } else {
      // Exact match required
      fRequiredMatch = requiredMatch;
    }
  }
  if (fRequiredMatch.empty()) {
    G4Exception
    ("G4PhysicalVolumesSearchScene::Matcher::Matcher",
     "modeling0013", JustWarning, "Required match is null");
  }
}

G4bool G4PhysicalVolumesSearchScene::Matcher::Match(const G4String& s)
// Match the string with the required match. The latter can be of the form
// "/regexp/", where regexp is a regular expression (see C++ regex),
// or a plain string, in which case there must be an exact match.
{
  G4bool found = false;
  if (fRequiredMatch.size()) {
    if (fRegexFlag) {  // Use extracted regex
      std::regex requiredRegex(fRequiredMatch);
      std::cmatch match;
      std::regex_search(s.c_str(), match, requiredRegex);
      if (match.size() > 0) found = true;
    } else {  // Require complete match
      if (s == fRequiredMatch) found = true;
    }
  }
  return found;
}
