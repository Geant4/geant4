#ifndef VOLUME_H
#define VOLUME_H 1

#include "SAXObject.hh"
#include "VolumeType.hh"

class volume : public SAXObject, public VolumeType
{
public:
  volume() {
  }
  virtual ~volume() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // VOLUME_H
