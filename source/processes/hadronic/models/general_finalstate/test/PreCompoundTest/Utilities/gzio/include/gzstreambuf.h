#ifndef gzstreambuf_h
#define gzstreambuf_h

#include <streambuf>
#include <zlib.h>

template <class charT, class traits = std::char_traits<charT> >
class basic_gzstreambuf : public std::basic_streambuf<charT,traits>
{
 public:
  typedef charT                            char_type;
  typedef traits                           traits_type;
  typedef typename traits_type::int_type   int_type;
  typedef typename traits_type::char_type  char_type;

  basic_gzstreambuf();
  ~basic_gzstreambuf();

  bool is_open() const;
  basic_gzstreambuf<charT,traits>* 
  open(const char* name, const std::ios_base::openmode open_mode, 
       const int level = GZ_DEFAULT_COMPRESSION);
  basic_gzstreambuf<charT, traits>* close();


 protected:
  std::streamsize xsputn(const char_type *s, std::streamsize n);
  int_type overflow(int_type c = traits_type::eof());
  int sync();
  int_type underflow();
  int_type pbackfail(int_type c);


 private:
  int flush_buffer();
  int fill_buffer();

  // prohibit copying and assignement
  basic_gzstreambuf(const basic_gzstreambuf & );
  basic_gzstreambuf & operator=(const basic_gzstreambuf & );
  
  //  static int PBSIZE() { return pbSize;}

 private:
  static const int           bufSize = 1048576;// 8192; //47+256;
  static const int           pbSize = 4;
  char_type                  buffer[bufSize];
  
  gzFile                     file;
  bool                       opened;
  std::ios_base::openmode    mode;
};

#include "gzstreambuf.icc"

typedef basic_gzstreambuf<char>     gzstreambuf;
typedef basic_gzstreambuf<wchar_t>  wgzstreambuf;


#endif // gzstreambuf_h
