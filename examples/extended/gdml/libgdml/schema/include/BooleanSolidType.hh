#ifndef BOOLEANSOLIDTYPE_H
#define BOOLEANSOLIDTYPE_H 1

#include "SAXObject.hh"

#include "SolidType.hh"
#include "ContentGroup.hh"
#include "ReferenceType.hh"
#include "positionType.hh"
#include "rotationType.hh"

class BooleanSolidType : public SolidType {
public:
  class first : public SAXObject, public ReferenceType {
  public:
    first() {
    }
    virtual ~first() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class second : public SAXObject, public ReferenceType {
  public:
    second() {
    }
    virtual ~second() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class position : public SAXObject, public positionType {
  public:
    position() {
    }
    virtual ~position() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class positionref : public SAXObject, public ReferenceType {
  public:
    positionref() {
    }
    virtual ~positionref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class rotation : public SAXObject, public rotationType {
  public:
    rotation() {
    }
    virtual ~rotation() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
  class rotationref : public SAXObject, public ReferenceType {
  public:
    rotationref() {
    }
    virtual ~rotationref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
public:
  BooleanSolidType() {
  }
  
  ~BooleanSolidType() {
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

#endif // BOOLEANSOLIDTYPE_H
