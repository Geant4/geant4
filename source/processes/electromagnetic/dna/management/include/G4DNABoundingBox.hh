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
// Author: HoangTRAN, 20/2/2019

#ifndef G4DNABoundingBox_hh
#define G4DNABoundingBox_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"
class G4DNABoundingBox
{
 public:
  G4DNABoundingBox()  = default;
  ~G4DNABoundingBox() = default;
  template <typename Iterator>
  G4DNABoundingBox(Iterator begin, Iterator end);
  G4DNABoundingBox(const std::initializer_list<G4double>& l);
  G4DNABoundingBox(G4DNABoundingBox&& rhs) noexcept;
  G4DNABoundingBox(const G4DNABoundingBox&) = default;
  G4DNABoundingBox& operator=(const G4DNABoundingBox&) = default;
  G4DNABoundingBox& operator=(G4DNABoundingBox&& hrs) noexcept;
  G4bool contains(const G4DNABoundingBox& other) const;
  G4bool contains(const G4ThreeVector& point) const;
  G4bool overlap(const G4DNABoundingBox& other, G4DNABoundingBox* out) const;
  G4bool overlap(const G4ThreeVector& query, const G4double& radius) const;
  G4bool contains(const G4ThreeVector& query, const G4ThreeVector& Point,
                  const G4double& Radius) const;
  G4bool contains(const G4ThreeVector& query, const G4double& radius) const;
  G4double Volume() const;
  std::array<G4DNABoundingBox, 8> partition() const;
  friend std::ostream& operator<<(std::ostream& s, const G4DNABoundingBox& rhs);
  G4bool operator==(const G4DNABoundingBox& rhs) const;
  G4bool operator!=(const G4DNABoundingBox& rhs) const;
  void resize(G4ThreeVector pics[8]);
  G4DNABoundingBox translate(const G4ThreeVector& trans) const;
  inline G4double halfSideLengthInX() const;
  inline G4double halfSideLengthInY() const;
  inline G4double halfSideLengthInZ() const;
  inline G4double Getxhi() const;
  inline G4double Getyhi() const;
  inline G4double Getzhi() const;
  inline G4double Getxlo() const;
  inline G4double Getylo() const;
  inline G4double Getzlo() const;
  inline G4ThreeVector middlePoint() const;

 private:
  G4double fxhi, fxlo, fyhi, fylo, fzhi, fzlo;
};

const G4DNABoundingBox initial = G4DNABoundingBox{
  -std::numeric_limits<G4double>::max(), std::numeric_limits<G4double>::max(),
  -std::numeric_limits<G4double>::max(), std::numeric_limits<G4double>::max(),
  -std::numeric_limits<G4double>::max(), std::numeric_limits<G4double>::max()
};
const G4DNABoundingBox invalid =
  G4DNABoundingBox{ std::numeric_limits<G4double>::quiet_NaN(),
                    std::numeric_limits<G4double>::quiet_NaN(),
                    std::numeric_limits<G4double>::quiet_NaN(),
                    std::numeric_limits<G4double>::quiet_NaN(),
                    std::numeric_limits<G4double>::quiet_NaN(),
                    std::numeric_limits<G4double>::quiet_NaN() };

inline G4ThreeVector G4DNABoundingBox::middlePoint() const
{
  G4ThreeVector mid((fxhi + fxlo) / 2.0, (fyhi + fylo) / 2.0,
                    (fzhi + fzlo) / 2.0);
  return mid;
}

inline G4double G4DNABoundingBox::halfSideLengthInX() const
{
  return std::abs(fxhi - fxlo) * 0.5;
}

inline G4double G4DNABoundingBox::halfSideLengthInY() const
{
  return std::abs(fyhi - fylo) * 0.5;
}

inline G4double G4DNABoundingBox::halfSideLengthInZ() const
{
  return std::abs(fzhi - fzlo) * 0.5;
}

template <typename Iterator>
G4DNABoundingBox::G4DNABoundingBox(Iterator begin, Iterator end)
  : G4DNABoundingBox(initial)
{
  for(; begin != end; ++begin)
  {
    const G4ThreeVector& point = (*begin);

    if(point.x() < fxlo)
    {
      fxlo = point.x();
    }
    if(point.x() > fxhi)
    {
      fxhi = point.x();
    }

    if(point.y() < fylo)
    {
      fylo = point.y();
    }
    if(point.y() > fyhi)
    {
      fyhi = point.y();
    }

    if(point.z() < fzlo)
    {
      fzlo = point.z();
    }
    if(point.z() > fzhi)
    {
      fzhi = point.z();
    }
  }
}

inline G4double G4DNABoundingBox::Getxhi() const { return fxhi; }
inline G4double G4DNABoundingBox::Getyhi() const { return fyhi; }
inline G4double G4DNABoundingBox::Getzhi() const { return fzhi; }
inline G4double G4DNABoundingBox::Getxlo() const { return fxlo; }
inline G4double G4DNABoundingBox::Getylo() const { return fylo; }
inline G4double G4DNABoundingBox::Getzlo() const { return fzlo; }
#endif
