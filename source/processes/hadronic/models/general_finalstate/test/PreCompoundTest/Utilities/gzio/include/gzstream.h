#ifndef gzstream_h
#define gzstream_h

#include <iostream>
#include "gzstreambuf.h"

//============= basic_gzstreambase ============================================

template <class charT, class traits = std::char_traits<charT> >
class basic_gzstreambase : virtual public std::basic_ios<charT,traits>
{
 protected:
  basic_gzstreambuf<charT,traits> buffer;
  
 public:
  basic_gzstreambase();
  basic_gzstreambase(const char* name, const std::ios_base::openmode open_mode,
                     const int level=Z_DEFAULT_COMPRESSION);
  ~basic_gzstreambase();
  
  void open(const char* name, const std::ios_base::openmode open_mode,
            const int level=Z_DEFAULT_COMPRESSION);
  void close();
  basic_gzstreambuf<charT,traits>* rdbuf();
  bool is_open() const;
};


//============= basic_igzstream ===============================================


template<class charT, class traits = std::char_traits<charT> >
class basic_igzstream : public basic_gzstreambase<charT,traits>,
                        public std::basic_istream<charT,traits>
{
 public:
  basic_igzstream() : std::basic_istream<charT,traits>(&buffer)
  {}
  
  basic_igzstream(const char* name)
    : basic_gzstreambase<charT,traits>(name, std::ios_base::in),
      std::basic_istream<charT, traits>(&buffer)
  {}
  
  ~basic_igzstream() {}
  
//   basic_gzstreambuf<charT, traits>* rdbuf()
//   {
//     return basic_gzstreambuf<CharT,traits>::rdbuf();
//   }

  void open(const char* name)
  {
    basic_gzstreambase<charT, traits>::open(name, std::ios_base::in);
  }

};

//============= basic_ogzstream ===============================================

template <class charT, class traits = std::char_traits<charT> >
class basic_ogzstream : public basic_gzstreambase<charT,traits>, 
                        public std::basic_ostream<charT,traits>
{
 public:
  basic_ogzstream() : std::basic_ostream<charT,traits>(&buffer)
  {}

  basic_ogzstream(const char* name, const int level=Z_DEFAULT_COMPRESSION)
    : basic_gzstreambase<charT, traits>(name, std::ios_base::out, level),
      std::basic_ostream<charT, traits>(&buffer)
  {}
  
  ~basic_ogzstream() {}
  
//   basic_gzstreambuf<charT,traits>* rdbud()
//   {
//     return basic_gzstreambase<charT,traits>::rdbuf();
//   }

  void open(const char* name, const int level=Z_DEFAULT_COMPRESSION)
  {
    basic_gzstreambase<charT,traits>::open(name, std::ios_base::out, level);
  }
};

#include "gzstream.icc"


typedef basic_igzstream<char>    igzstream;
typedef basic_ogzstream<char>    ogzstream;
typedef basic_igzstream<wchar_t> wigzstream;
typedef basic_ogzstream<wchar_t> wogzstream;


#endif // gzstream_h
