#ifndef VOLUMETYPE_H
#define VOLUMETYPE_H 1

#include "IdentifiableVolumeType.hh"
#include "ReferenceType.hh"
#include "SinglePlacementType.hh"
#include "ContentGroup.hh"

class VolumeType : public IdentifiableVolumeType {
public:
    
  class materialref : public SAXObject, public ReferenceType {
  public:
    materialref() {
    }
    virtual ~materialref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class solidref : public SAXObject, public ReferenceType {
  public:
    solidref() {
    }
    virtual ~solidref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
public:
  VolumeType() {
  }
  ~VolumeType() {
  }
  
  const ContentSequence* get_content() const {
    return &m_sequence;
  }

  void add_content( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = { tag, so };
    m_sequence.add_content( ci );
  }
private:
  ContentSequence m_sequence;
};

#endif // VOLUMETYPE_H
