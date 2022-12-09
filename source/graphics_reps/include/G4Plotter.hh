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
// Guy Barrand 25th September 2021
//

#ifndef G4PLOTTER_HH
#define G4PLOTTER_HH

#include "G4String.hh"

#include <utility>
#include <vector>

namespace tools {namespace histo {class h1d;}}
namespace tools {namespace histo {class h2d;}}

class G4Plotter {
public:  
  using RegionStyle     = std::pair<unsigned int,G4String>;
  using Parameter       = std::pair<G4String,G4String>;
  using RegionParameter = std::pair<unsigned int,Parameter>;
  using Region_h1d      = std::pair<unsigned int,tools::histo::h1d*>;
  using Region_h2d      = std::pair<unsigned int,tools::histo::h2d*>;
  using Region_h1       = std::pair<unsigned int,int>;
  using Region_h2       = std::pair<unsigned int,int>;

  G4Plotter();
  virtual ~G4Plotter() = default;
  G4Plotter(const G4Plotter&);
  G4Plotter& operator = (const G4Plotter&);

  void SetLayout(unsigned int colums,unsigned int rows);
  void AddStyle(const G4String& style);
  void AddRegionStyle(unsigned int region,const G4String& style);
  void AddRegionParameter(unsigned int region,const G4String& parameter,const G4String& value);
  void AddRegionHistogram(unsigned int region,tools::histo::h1d* histo);
  void AddRegionHistogram(unsigned int region,tools::histo::h2d* histo);
  void AddRegionH1(unsigned int region,int id);
  void AddRegionH2(unsigned int region,int id);
  void Reset();
  void Clear();
  void ClearRegion(unsigned int region);

  unsigned int GetColumns() const {return fColumns;}
  unsigned int GetRows() const {return fRows;}
  const std::vector<G4String>& GetStyles() const {return fStyles;}
  const std::vector<RegionStyle>& GetRegionStyles() const {return fRegionStyles;}
  const std::vector<RegionParameter>& GetRegionParameters() const {return fRegionParameters;}
  const std::vector<Region_h1d>& GetRegionH1Ds() const {return fRegion_h1ds;}
  const std::vector<Region_h2d>& GetRegionH2Ds() const {return fRegion_h2ds;}

  const std::vector<Region_h1>& GetRegionH1s() const {return fRegion_h1s;}
  const std::vector<Region_h2>& GetRegionH2s() const {return fRegion_h2s;}

private:
  unsigned int fColumns{1};
  unsigned int fRows{1};
  std::vector<G4String> fStyles;
  std::vector<RegionStyle> fRegionStyles;
  std::vector<RegionParameter> fRegionParameters;
  std::vector<Region_h1d> fRegion_h1ds;
  std::vector<Region_h2d> fRegion_h2ds;
  std::vector<Region_h1>  fRegion_h1s;
  std::vector<Region_h2>  fRegion_h2s;
};

#endif
