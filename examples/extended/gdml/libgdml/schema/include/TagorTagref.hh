#ifndef TAGORTAGREF_H
#define TAGORTAGREF_H 1

#include "ContentGroup.hh"

class TagorTagref
{
public:
  TagorTagref() {
  }
  ~TagorTagref() {
  }

  const SAXObject* get_content() const {
    return &m_choice;
  }

  void set_content( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = {tag,so};
    m_choice.set_content( ci );
  }

private:
  // choice content model: ( Tag | ref )
  ContentChoice m_choice;
};

#endif // TAGORTAGREF_H
