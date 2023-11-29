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
#include "G4DNABoundingBox.hh"
#include <array>
#include "G4UnitsTable.hh"
using std::array;
using std::initializer_list;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNABoundingBox::G4DNABoundingBox(G4DNABoundingBox&& rhs) noexcept
{
  fxhi = rhs.fxhi;
  fxlo = rhs.fxlo;
  fyhi = rhs.fyhi;
  fylo = rhs.fylo;
  fzhi = rhs.fzhi;
  fzlo = rhs.fzlo;
}
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNABoundingBox::G4DNABoundingBox(const initializer_list<G4double>& l)
{
  fxhi = 0.;
  fxlo = 0.;
  fyhi = 0.;
  fylo = 0.;
  fzhi = 0.;
  fzlo = 0.;
  std::copy(l.begin(), l.end(), &fxhi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNABoundingBox& G4DNABoundingBox::operator=(G4DNABoundingBox&& rhs) noexcept
{
  fxhi = rhs.fxhi;
  fxlo = rhs.fxlo;
  fyhi = rhs.fyhi;
  fylo = rhs.fylo;
  fzhi = rhs.fzhi;
  fzlo = rhs.fzlo;
  return *this;
}

G4double G4DNABoundingBox::Volume() const
{
  auto LengthY = this->Getyhi() - this->Getylo();
  auto LengthX = this->Getxhi() - this->Getxlo();
  auto LengthZ = this->Getzhi() - this->Getzlo();
  G4double V   = LengthY * LengthX * LengthZ;
  return V;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNABoundingBox::resize(G4ThreeVector pics[8])
{
  for(size_t i = 0; i < 8; i++)
  {
    const G4ThreeVector& point = pics[i];
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

//------------------------------------------------------------------------------

G4DNABoundingBox G4DNABoundingBox::translate(const G4ThreeVector& trans) const
{
  G4DNABoundingBox output;

  output.fxhi = fxhi + trans.x();
  output.fxlo = fxlo + trans.x();

  output.fyhi = fyhi + trans.y();
  output.fylo = fylo + trans.x();

  output.fzhi = fzhi + trans.z();
  output.fzlo = fzlo + trans.z();

  return output;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::contains(const G4DNABoundingBox& other) const
{
  return fxlo <= other.fxlo && fxhi >= other.fxhi && fylo <= other.fylo &&
         fyhi >= other.fyhi && fzlo <= other.fzlo && fzhi >= other.fzhi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::contains(const G4ThreeVector& point) const
{
  return fxlo <= point.x() && fxhi >= point.x() && fylo <= point.y() &&
         fyhi >= point.y() && fzlo <= point.z() && fzhi >= point.z();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::overlap(const G4DNABoundingBox& other,
                                 G4DNABoundingBox* out) const
{
  if(contains(other))
  {
    *out = other;
    return true;
  }
  else if(other.contains(*this))
  {
    *out = *this;
    return true;
  }

  // Check if there is no intersection
  if(fxhi < other.fxlo || fxlo > other.fxhi || fyhi < other.fylo ||
     fylo > other.fyhi || fzhi < other.fzlo || fzlo > other.fzhi)
  {
    *out = invalid;
    return false;
  }

  G4double upperX = std::min(fxhi, other.fxhi);
  G4double upperY = std::min(fyhi, other.fyhi);
  G4double upperZ = std::min(fzhi, other.fzhi);

  G4double lowerX = std::max(fxlo, other.fxlo);
  G4double lowerY = std::max(fylo, other.fylo);
  G4double lowerZ = std::max(fzlo, other.fzlo);

  *out = G4DNABoundingBox{ upperX, lowerX, upperY, lowerY, upperZ, lowerZ };
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::overlap(const G4ThreeVector& query,
                                 const G4double& Radius) const
{
  G4double x = query.x() - this->middlePoint().x();
  G4double y = query.y() - this->middlePoint().y();
  G4double z = query.z() - this->middlePoint().z();

  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);

  if((x > (Radius + this->halfSideLengthInX())) ||
     (y > (Radius + this->halfSideLengthInY())) ||
     (z > (Radius + this->halfSideLengthInZ())))
  {
    return false;  //
  }
  G4int num_less_extent = (x < this->halfSideLengthInX()) +
                          (y < this->halfSideLengthInY()) +
                          (z < this->halfSideLengthInZ());

  if(num_less_extent > 1)
  {
    return true;
  }

  x = std::max(x - this->halfSideLengthInX(), 0.0);
  y = std::max(y - this->halfSideLengthInY(), 0.0);
  z = std::max(z - this->halfSideLengthInZ(), 0.0);

  G4double norm = std::sqrt(x * x + y * y + z * z);

  return (norm < Radius);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::contains(const G4ThreeVector& query,
                                  const G4double& Radius) const
{
  G4ThreeVector temp = query - this->middlePoint();
  G4double x         = std::abs(temp.x()) + this->halfSideLengthInX();
  G4double y         = std::abs(temp.y()) + this->halfSideLengthInY();
  G4double z         = std::abs(temp.z()) + this->halfSideLengthInZ();
  G4double norm      = std::sqrt(x * x + y * y + z * z);
  return (norm < Radius);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::contains(const G4ThreeVector& query,
                                  const G4ThreeVector& Point,
                                  const G4double& Radius) const
{
  return (((query - Point).mag()) < Radius);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

array<G4DNABoundingBox, 8> G4DNABoundingBox::partition() const
{
  G4double xmid = (fxhi + fxlo) / 2.;
  G4double ymid = (fyhi + fylo) / 2.;
  G4double zmid = (fzhi + fzlo) / 2.;

  std::array<G4DNABoundingBox, 8> ret{ {
    G4DNABoundingBox{ xmid, fxlo, ymid, fylo, zmid,
                      fzlo },  // bottom left front
    G4DNABoundingBox{ fxhi, xmid, ymid, fylo, zmid,
                      fzlo },  // bottom right front
    G4DNABoundingBox{ xmid, fxlo, fyhi, ymid, zmid, fzlo },  // bottom left back
    G4DNABoundingBox{ fxhi, xmid, fyhi, ymid, zmid,
                      fzlo },  // bottom right back
    G4DNABoundingBox{ xmid, fxlo, ymid, fylo, fzhi, zmid },  // top left front
    G4DNABoundingBox{ fxhi, xmid, ymid, fylo, fzhi, zmid },  // top right front
    G4DNABoundingBox{ xmid, fxlo, fyhi, ymid, fzhi, zmid },  // top left back
    G4DNABoundingBox{ fxhi, xmid, fyhi, ymid, fzhi, zmid }   // top right back
  } };
  return ret;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4DNABoundingBox::operator==(const G4DNABoundingBox& rhs) const
{
  return (fxhi == rhs.fxhi && fxlo == rhs.fxlo && fyhi == rhs.fyhi &&
          fylo == rhs.fylo && fzhi == rhs.fzhi && fzlo == rhs.fzlo) ||
         (std::isnan(fxhi) && std::isnan(rhs.fxhi) && std::isnan(fxlo) &&
          std::isnan(rhs.fxlo) && std::isnan(fyhi) && std::isnan(rhs.fyhi) &&
          std::isnan(fylo) && std::isnan(rhs.fylo) && std::isnan(fzhi) &&
          std::isnan(rhs.fzhi) && std::isnan(fzlo) && std::isnan(rhs.fzlo));
}

G4bool G4DNABoundingBox::operator!=(const G4DNABoundingBox& rhs) const
{
  return !operator==(rhs);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream& operator<<(std::ostream& stream, const G4DNABoundingBox& rhs)
{
  stream << "{" << G4BestUnit(rhs.fxhi, "Length") << ", "
         << G4BestUnit(rhs.fxlo, "Length") << ", "
         << G4BestUnit(rhs.fyhi, "Length") << ", "
         << G4BestUnit(rhs.fylo, "Length") << ", "
         << G4BestUnit(rhs.fzhi, "Length") << ", "
         << G4BestUnit(rhs.fzlo, "Length") << ", "
         << "}";
  return stream;
}
