#ifndef SINGLEPLACEMENTTYPE_H
#define SINGLEPLACEMENTTYPE_H 1

#include "ContentGroup.hh"
#include "ReferenceType.hh"

class SinglePlacementType {
public:
    
  class volumeref : public SAXObject, public ReferenceType {
  public:
    volumeref() {
    }
    virtual ~volumeref() {
    }
    virtual SAXObject::Type type() {
      return SAXObject::element;
    }
  };
  
//  We need to resolve the problem of conlifting non-terminals in grammar
//  This is typical use-case here, the following two elements clash with
//  BooleanSolidType elements
//   class positionref : public SAXObject, public ReferenceType {
//   public:
//     positionref() {
//     }
//     virtual ~positionref() {
//     }
//     virtual SAXObject::Type type() {
//       return SAXObject::element;
//     }
//   };
//   
//   class rotationref : public SAXObject, public ReferenceType {
//   public:
//     rotationref() {
//     }
//     virtual ~rotationref() {
//     }
//     virtual SAXObject::Type type() {
//       return SAXObject::element;
//     }
//   };
  
public:
  SinglePlacementType() {
  }
  ~SinglePlacementType() {
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

#endif // SINGLEPLACEMENTTYPE_H
