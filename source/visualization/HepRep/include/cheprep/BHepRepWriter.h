// Copyright FreeHEP, 2005.
#ifndef CHEPREP_BHEPREPWRITER_H
#define CHEPREP_BHEPREPWRITER_H

#include "cheprep/config.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "cheprep/AbstractXMLWriter.h"

/**
 * @author Mark Donszelmann
 * @version $Id: BHepRepWriter.h 66870 2013-01-14 23:38:59Z adotti $
 */
namespace cheprep {

    class BHepRepWriter : public AbstractXMLWriter {

        public:

            BHepRepWriter(std::ostream& os);
            virtual ~BHepRepWriter();

            void close();
            void openDoc(std::string version = "BinaryHepRep/1.0", std::string encoding = "UTF-8", bool standalone = false);
            void closeDoc(bool force = false);
            void openTag(std::string name);
            void closeTag();
            void printTag(std::string name);
            void setAttribute(std::string name, char* value);
            void setAttribute(std::string name, std::string value);
            void setAttribute(std::string name, std::vector<double> value);
            void setAttribute(std::string name, int64 value);
            void setAttribute(std::string name, int value);
            void setAttribute(std::string name, bool value);
            void setAttribute(std::string name, double value);

            //
            // Can be removed when we can properly inherit those (since names are equal to overloaded ones).
            // 
            void openTag(std::string ns, std::string name) {
                openTag(ns == defaultNameSpace ? name : ns.append(":").append(name));
            }
            void printTag(std::string ns, std::string name) {
                printTag(ns == defaultNameSpace ? name : ns.append(":").append(name));
            }
            void setAttribute(std::string ns, std::string name, std::string value) {
                setAttribute(ns.append(":").append(name), value);
            }
            void setAttribute(std::string ns, std::string name, double value) {
                setAttribute(ns.append(":").append(name), value);
            }

        private:
            static const unsigned char WBXML_VERSION    = 0x03;
            static const unsigned char UNKNOWN_PID      = 0x01;
            static const unsigned char UTF8             = 0x6a;
        
            // standard tags
            static const unsigned char SWITCH_PAGE      = 0x00;
            static const unsigned char END              = 0x01;
            static const unsigned char ENTITY           = 0x02;
            static const unsigned char STR_I            = 0x03;
            static const unsigned char LITERAL          = 0x04;
            
            static const unsigned char CONTENT          = 0x40;
            static const unsigned char EXT_I_0          = 0x40;
            static const unsigned char EXT_I_1          = 0x41;
            static const unsigned char EXT_I_2          = 0x42;
            static const unsigned char PI               = 0x43;
            static const unsigned char LITERAL_C        = 0x44;
            
            static const unsigned char ATTRIBUTE        = 0x80;
            static const unsigned char EXT_T_0          = 0x80;
            static const unsigned char EXT_T_1          = 0x81;
            static const unsigned char EXT_T_2          = 0x82;
            static const unsigned char STR_T            = 0x83;
            static const unsigned char LITERAL_A        = 0x84;
            
            static const unsigned char EXT_0            = 0xC0;
            static const unsigned char EXT_1            = 0xC1;
            static const unsigned char EXT_2            = 0xC2;
            static const unsigned char OPAQUE           = 0xC3;
            static const unsigned char LITERAL_AC       = 0xC4;
            
            // our own extensions            
            static const unsigned char STR_D            = EXT_I_0;
            static const unsigned char STR_R            = EXT_T_0;
            
            // class definitions
            static std::map<std::string, unsigned char> tags; 
            static std::map<std::string, unsigned char> attributes;           
            static std::map<std::string, unsigned char> values;

            // outputstream variables
            std::ostream& os;
            bool singlePrecision;    
            bool isBigEndian;
            
            // document variables            
            std::map<std::string, unsigned int> stringValues;
            
            // tag variables
            std::map<std::string, std::string> stringAttributes;
            std::map<std::string, std::vector<double> > colorAttributes;
            std::map<std::string, int64> longAttributes;
            std::map<std::string, int> intAttributes;
            std::map<std::string, bool> booleanAttributes;
            std::map<std::string, double> doubleAttributes;
            
            // point array
            std::vector<double> points;
            
            // methods
            void writeTag(std::string name, bool content = false);
            void writePoints();
            void writeStringDefine(std::string s);                     
            void writeMultiByteInt(unsigned int ui);                     
            void writeReal(double ui);
            void writeLong(int64 i);
            void writeInt(int i);
            void writeByte(unsigned char b);
            void writeString(std::string s);
    };
 
} // cheprep

#endif // CHEPREP_BHEPREPWRITER_H
