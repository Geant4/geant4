#ifndef ASSEMBLY_H
#define ASSEMBLY_H 1

#include "SAXObject.hh"
#include "AssemblyVolumeType.hh"

class assembly : public SAXObject, public AssemblyVolumeType
{
public:
  assembly() {
  }
  virtual ~assembly() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // ASSEMBLY_H
