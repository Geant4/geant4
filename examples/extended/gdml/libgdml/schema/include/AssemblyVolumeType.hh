#ifndef ASSEMBLYVOLUMETYPE_H
#define ASSEMBLYVOLUMETYPE_H 1

#include "IdentifiableVolumeType.hh"
#include "SinglePlacementType.hh"
#include "ContentGroup.hh"

class AssemblyVolumeType : public IdentifiableVolumeType {
public:
  AssemblyVolumeType() {
  }
  ~AssemblyVolumeType() {
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

#endif // ASSEMBLYVOLUMETYPE_H
