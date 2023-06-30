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

#ifndef G4VVTKPIPELINE_HH
#define G4VVTKPIPELINE_HH

#include "G4Polyhedron.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4VisAttributes.hh"
#include "G4VtkVisContext.hh"

#include "vtkRenderer.h"

namespace std
{
inline void hash_combine(std::size_t) {}

template<typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  std::hash_combine(seed, rest...);
}

template<>
struct hash<G4String>
{
    std::size_t operator()(const G4String& string) const
    {
      using std::hash;
      using std::size_t;

      std::size_t h = 0;

      for (char const& c : string) {
        std::hash_combine(h, c);
      }

      return h;
    }
};

template<>
struct hash<G4VisAttributes>
{
    std::size_t operator()(const G4VisAttributes& va) const
    {
      using std::hash;
      using std::size_t;

      std::size_t h = 0;

      std::hash_combine(h, va.IsVisible());
      std::hash_combine(h, va.IsDaughtersInvisible());
      std::hash_combine(h, va.GetColour().GetRed());
      std::hash_combine(h, va.GetColour().GetGreen());
      std::hash_combine(h, va.GetColour().GetBlue());
      std::hash_combine(h, va.GetColour().GetAlpha());
      std::hash_combine(h, static_cast<int>(va.GetLineStyle()));

      return h;
    }
};

template<>
struct hash<G4Polyhedron>
{
    std::size_t operator()(const G4Polyhedron& ph) const
    {
      using std::hash;
      using std::size_t;

      G4bool notLastFace;
      G4Point3D vertex[4];
      G4int edgeFlag[4];
      G4Normal3D normals[4];
      G4int nEdges;

      std::size_t h = 0;

      do {
        notLastFace = ph.GetNextFacet(nEdges, vertex, edgeFlag, normals);

        for (int i = 0; i < nEdges; i++) {
          std::size_t hx = std::hash<double>()(vertex[i].x());
          std::size_t hy = std::hash<double>()(vertex[i].y());
          std::size_t hz = std::hash<double>()(vertex[i].z());
          std::hash_combine(h, hx);
          std::hash_combine(h, hy);
          std::hash_combine(h, hz);
        }
      } while (notLastFace);

      return h;
    }
};
}  // namespace std

class G4VVtkPipeline
{
  public:
    G4VVtkPipeline() : name("none"), type("G4VVtkPipeline"), disableParent(false), renderer(nullptr)
    {}
    G4VVtkPipeline(G4String nameIn, G4String typeIn, const G4VtkVisContext& vcIn,
                   G4bool disableParentIn = false,
                   vtkSmartPointer<vtkRenderer> rendererIn = nullptr)
      : name(nameIn), type(typeIn), disableParent(disableParentIn), renderer(rendererIn), vc(vcIn)
    {}
    virtual ~G4VVtkPipeline()
    {
      for (auto i : childPipelines)
        delete i;
    }

    virtual void Enable() = 0;
    virtual void Disable() = 0;

    virtual void Print()
    {
      G4cout << "children (";
      for (auto i : childPipelines)
        G4cout << i->GetName() << "-" << i->GetTypeName() << ",";
      G4cout << ")" << G4endl;
    };
    virtual void Modified()
    {
      for (auto i : childPipelines)
        i->Modified();
    }
    virtual void Clear()
    {
      for (auto i : childPipelines)
        i->Clear();
    }

    void SetDisableParent(G4bool disableParentIn) { disableParent = disableParentIn; }
    bool GetDisableParent() { return disableParent; }

    void SetName(G4String nameIn) { name = nameIn; }
    const G4String GetName() { return name; }

    void SetTypeName(G4String typeNameIn) { type = typeNameIn; }
    G4String GetTypeName() { return type; }

    G4VtkVisContext& GetVtkVisContext() { return vc; }

    void AddChildPipeline(G4VVtkPipeline* child)
    {
      childPipelines.push_back(child);
      if (child->GetDisableParent()) {
        Disable();
      }
    }
    G4VVtkPipeline* GetChildPipeline(G4int iPipeline) { return childPipelines[iPipeline]; }
    G4int GetNumberOfChildPipelines() { return (G4int)childPipelines.size(); }
    std::vector<G4VVtkPipeline*> GetChildPipelines() { return childPipelines; }
    void ClearChildPipeline() { childPipelines.clear(); }

  protected:
    G4String name;
    G4String type;
    G4bool disableParent;
    std::vector<G4VVtkPipeline*> childPipelines;
    vtkSmartPointer<vtkRenderer> renderer;
    G4VtkVisContext vc;
};

#endif  // G4VVTKPIPELINE_HH
