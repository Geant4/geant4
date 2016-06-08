#ifndef rwstlhash_h
#define rwstlhash_h

class G4ParticleDefinition;
class _WidgetRec;

inline unsigned HashDefault( const G4ParticleDefinition* const & p)
{
  return (long) p;
}

inline unsigned HashDefault( _WidgetRec* const & p)
{
  return (long) p;
}

inline unsigned HashDefault( const int& p) 
{
  return p;
}

template < class K > 
unsigned HashDefault( const K& p)
{
  return p.stlhash();
}


#endif
