#ifndef MATERIALPROPERTIESGROUP_H
#define MATERIALPROPERTIESGROUP_H 1

#include "TagorTagref.hh"


// Morphism pattern used:
// element group -> struct
// choice        -> ContentChoice      (TagorTagref)
// choice?       -> ContentChoice* ptr (TagorTagref* ptr)
struct  MaterialPropertiesGroup
{
  // optional RL or RLref
  TagorTagref* m_RLorRLref;
  // optional AL or ALref
  TagorTagref* m_ALorALref;
  // optional T or Tref
  TagorTagref* m_TorTref;
  // optional P or Pref
  TagorTagref* m_PorPref;

  MaterialPropertiesGroup()
  : m_RLorRLref(0), m_ALorALref(0), m_TorTref(0), m_PorPref(0) {
  }
};

#endif // MATERIALPROPERTIESGROUP_H
