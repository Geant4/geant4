// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFLATEOUTPUTSTREAMBUF_H
#define CHEPREP_DEFLATEOUTPUTSTREAMBUF_H

#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

#ifndef CHEPREP_NO_ZLIB
#include <zlib.h>
#endif // CHEPREP_NO_ZLIB

/**
 * @author Mark Donszelmann
 * @version $Id: DeflateOutputStreamBuffer.h 66870 2013-01-14 23:38:59Z adotti $
 */
namespace cheprep {

    class DeflateOutputStreamBuffer : public std::streambuf {

        public:
        
            DeflateOutputStreamBuffer(std::streambuf *buffer);
        
            void init(bool compress);
            void finish();
        
            virtual ~DeflateOutputStreamBuffer();
    
        
        protected:
            int overflow(int c = EOF);

#ifndef CHEPREP_NO_ZLIB
            bool flushOut();
#endif // CHEPREP_NO_ZLIB

            inline void putUI(unsigned int ui) {
                unsigned char* ucp = reinterpret_cast<unsigned char *>(&ui);
                unsigned int i   = (static_cast<unsigned int>(ucp[ 3 ]) << 24) +  
                                   (static_cast<unsigned int>(ucp[ 2 ]) << 16) + 
                                   (static_cast<unsigned int>(ucp[ 1 ]) << 8 ) + 
                                   (static_cast<unsigned int>(ucp[ 0 ]));                   
                buffer->sputn(reinterpret_cast<char *>(&i), sizeof(unsigned int));
            }

            inline void putUS(unsigned short us) {
                unsigned char* ucp = reinterpret_cast<unsigned char *>(&us);
                unsigned short s   = (static_cast<unsigned short>(ucp[ 1 ]) << 8 ) + 
                                     (static_cast<unsigned short>(ucp[ 0 ]));                   
                buffer->sputn(reinterpret_cast<char *>(&s), sizeof(unsigned short));
            }

            inline void putUB(unsigned char ub) {
                buffer->sputc(ub);
            }

            inline void putS(const std::string s) {
                buffer->sputn(s.c_str(), s.length());                
            }
            
            inline std::streampos pos() {
                std::ostream os(buffer);
                return os.tellp();
            }
            
            inline unsigned int getSize() {
                return size;
            }
            
            inline unsigned int getCRC() {
                return crc;
            }
            
        private:
            static unsigned long crctable[256];                
            std::streambuf *buffer;
            
            unsigned int crc;
            unsigned int size;

#ifndef CHEPREP_NO_ZLIB            
            static const unsigned int inSize = 1000;
            static const unsigned int outSize = 1000;
            z_stream zStream;
            bool zStreamOpen;

            std::vector<char> in;
            std::vector<char> out;
#endif // CHEPREP_NO_ZLIB
};


} // cheprep



#endif // CHEPREP_DEFLATEOUTPUTSTREAMBUF_H 
