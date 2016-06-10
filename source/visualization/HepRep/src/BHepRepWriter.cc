// Copyright FreeHEP, 2005.

#include "cheprep/config.h"

#include <iostream>
#include <algorithm>

#include "cheprep/BHepRepWriter.h"
                        
/**
 * @author Mark Donszelmann
 * @version $Id: BHepRepWriter.cc 66870 2013-01-14 23:38:59Z adotti $
 */
namespace cheprep {

    // class variables
    std::map<std::string, unsigned char> BHepRepWriter::tags; 
    std::map<std::string, unsigned char> BHepRepWriter::attributes;           
    std::map<std::string, unsigned char> BHepRepWriter::values;

    BHepRepWriter::BHepRepWriter( std::ostream& ostrm) 
            : AbstractXMLWriter("heprep"), 
            os(ostrm), 
            singlePrecision(true) {        

        // resolve endiannes        
        union { long l; char c[sizeof (long)]; } u;
        u.l = 1;
        isBigEndian = (u.c[sizeof (long) - 1] == 1);
               
//        std::cout << "Host is " << (isBigEndian ? "Big-Endian" : "Little-Endian") << "." << std::endl;
       
        if (tags.size() <= 0) { 
            // tags
            tags["heprep"]              = 0x05;
            tags["attdef"]              = 0x06;
            tags["attvalue"]            = 0x07;
            tags["instance"]            = 0x08;
            tags["treeid"]              = 0x09;
            tags["action"]              = 0x0a;
            tags["instancetree"]        = 0x0b;
            tags["type"]                = 0x0c;
            tags["typetree"]            = 0x0d;
            tags["layer"]               = 0x0e;
            tags["point"]               = 0x0f;
        }

        if (attributes.size() <= 0) {
            // attribute names
            attributes["version"]               = 0x05;
            attributes["xmlns"]                 = 0x06;
            attributes["xmlns:xsi"]             = 0x07;
            attributes["xsi:schemaLocation"]    = 0x08;
    
            attributes["valueString"]           = 0x10;
            attributes["valueColor"]            = 0x11;
            attributes["valueLong"]             = 0x12;
            attributes["valueInt"]              = 0x13;
            attributes["valueBoolean"]          = 0x14;
            attributes["valueDouble"]           = 0x15;
    
            attributes["name"]                  = 0x20;
            attributes["type"]                  = 0x22;
            attributes["showlabel"]             = 0x23;
            attributes["desc"]                  = 0x24;
            attributes["category"]              = 0x25;
            attributes["extra"]                 = 0x26;
            attributes["x"]                     = 0x27;
            attributes["y"]                     = 0x28;
            attributes["z"]                     = 0x29;
            attributes["qualifier"]             = 0x2a;
            attributes["expression"]            = 0x2b;
            attributes["typetreename"]          = 0x2c;
            attributes["typetreeversion"]       = 0x2d;
            attributes["order"]                 = 0x2e;
            
            // for PI
            attributes["eof"]                   = 0x7f;
        }
        
        if (values.size() <= 0) { 
            // attribute values
            values["drawas"]                     = 0x85;
            values["drawasoptions"]              = 0x86;
            values["visibility"]                 = 0x87;
    
            values["label"]                      = 0x88;
    
            values["fontname"]                   = 0x89;
            values["fontstyle"]                  = 0x8a;
            values["fontsize"]                   = 0x8b;
            values["fontcolor"]                  = 0x8c;
            values["fonthasframe"]               = 0x8d;
            values["fontframecolor"]             = 0x8e;
            values["fontframewidth"]             = 0x8f;
            values["fonthasbanner"]              = 0x90;
            values["fontbannercolor"]            = 0x91;
    
            values["color"]                      = 0x92;
            values["framecolor"]                 = 0x93;
            values["layer"]                      = 0x94;
            values["markname"]                   = 0x95;
            values["marksize"]                   = 0x96;
            values["marksizemultiplier"]         = 0x97;
            values["marktype"]                   = 0x98;
            values["hasframe"]                   = 0x99;
            values["framecolor"]                 = 0x9a;
            values["framewidth"]                 = 0x9b;
    
            values["linestyle"]                  = 0x9c;
            values["linewidth"]                  = 0x9d;
            values["linewidthmultiplier"]        = 0x9e;
            values["linehasarrow"]               = 0x9f;
            
            values["fillcolor"]                  = 0xa0;
            values["filltype"]                   = 0xa1;
            values["fill"]                       = 0xa2;
    
            values["radius"]                     = 0xa3;
            values["phi"]                        = 0xa4;
            values["theta"]                      = 0xa5;
            values["omega"]                      = 0xa6;
            values["radius1"]                    = 0xa7;
            values["radius2"]                    = 0xa8;
            values["radius3"]                    = 0xa9;
            values["curvature"]                  = 0xaa;
            values["flylength"]                  = 0xab;
            values["faces"]                      = 0xac;
    
            values["text"]                       = 0xad;
            values["hpos"]                       = 0xae;
            values["vpos"]                       = 0xaf;
            values["halign"]                     = 0xb0;
            values["valign"]                     = 0xb1;
    
            values["ispickable"]                 = 0xb2;
            values["showparentvalues"]           = 0xb3;
            values["pickparent"]                 = 0xb4;
            
            // attvalue values
            values["false"]         = 0xd0;
            values["true"]          = 0xd1;
            
            values["point"]         = 0xd2;
            values["line"]          = 0xd3;
            values["helix"]         = 0xd4;
            values["polygon"]       = 0xd5;
            values["circle"]        = 0xd6;
            values["curve"]         = 0xd7;
            values["ellipse"]       = 0xd8;
            values["ellipsoid"]     = 0xd9;
            values["prism"]         = 0xda;
            values["cylinder"]      = 0xdb;
            values["ellipseprism"]  = 0xdc;
            values["text"]          = 0xdd;
    
            values["nonzero"]       = 0xde;
            values["evenodd"]       = 0xdf;
            
            values["circle"]        = 0xe0;
            values["box"]           = 0xe1;
            values["uptriangle"]    = 0xe2;
            values["dntriangle"]    = 0xe3;
            values["diamond"]       = 0xe4;
            values["cross"]         = 0xe5;
            values["star"]          = 0xe6;
            values["plus"]          = 0xe7;
            values["hline"]         = 0xe8;
            values["vline"]         = 0xe9;
    
            values["solid"]         = 0xea;
            values["dotted"]        = 0xeb;
            values["dashed"]        = 0xec;
            values["dotdash"]       = 0xed;
            
            values["none"]          = 0xee;
            values["start"]         = 0xef;
            values["end"]           = 0xf0;
            values["both"]          = 0xf1;
    
            values["serif"]         = 0xf2;
            values["sansserif"]     = 0xf3;
            values["monotype"]      = 0xf4;
            values["symbol"]        = 0xf5;
    
            values["plain"]         = 0xf6;
            values["bold"]          = 0xf7;
            values["italic"]        = 0xf8;
    
            values["top"]           = 0xf9;
            values["baseline"]      = 0xfa;
            values["center"]        = 0xfb;
            values["bottom"]        = 0xfc;
    
            values["left"]          = 0xfd;
            values["right"]         = 0xfe;
    
            values["default"]       = 0xff;
        }        
    }
    
    BHepRepWriter::~BHepRepWriter() {
    }

    void BHepRepWriter::close() {
    }
    
    void BHepRepWriter::openDoc(std::string version, std::string /* encoding */, bool /* standalone */) {
        stringValues.clear();
        
        // header
        writeByte(WBXML_VERSION);
        writeMultiByteInt(UNKNOWN_PID);
        writeMultiByteInt(UTF8);        
        
        version = "BinaryHepRep/1.0"; 
       
        // string table
        writeMultiByteInt(version.length()+1);
        
        // BHepRep Header (as part of the string table)
        writeString(version);
        
    }
    
    void BHepRepWriter::closeDoc(bool /* force */) {
        writeByte(PI);
        writeByte(attributes["eof"]);
        writeByte(END);
    }

    void BHepRepWriter::openTag(std::string name) {
        writeTag(name, true);
    }
    
    void BHepRepWriter::closeTag() {
        writePoints();
        writeByte(END);
    }
    
    void BHepRepWriter::printTag(std::string name) {
        writeTag(name);
    }
    
    void BHepRepWriter::writeTag(std::string tagName, bool hasContent) {
        std::string s = tagName;
        std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
        
        // find tag
        if (tags.count(s) <= 0) {
            std::cerr << "Cannot find tag '" << s << "' in tags table." << std::endl;
            return;
        }
        
        // write tag
        bool isPoint = (s == "point");
        bool hasAttributes = (stringAttributes.size() > 0) || (doubleAttributes.size() > (unsigned int)(isPoint ? 3 : 0));
        
        if (!hasAttributes && isPoint) {
            // store the point for the future
            points.push_back(doubleAttributes["x"]);
            points.push_back(doubleAttributes["y"]);
            points.push_back(doubleAttributes["z"]);
            return;
        }

        writePoints();
        writeByte(tags[s] | ((hasContent || isPoint) ? CONTENT : 0x00) | (hasAttributes ? ATTRIBUTE : 0x00));        
            
        // write attributes
        if (hasAttributes) {
            // write string attributes
	        for (std::map<std::string,std::string>::iterator i = stringAttributes.begin(); i != stringAttributes.end(); i++) {
    		    std::string name = i->first;
    		    std::string value = i->second;

                // write ATTRSTART
                writeByte(attributes[name]);
                std::string v = value;
                std::transform(v.begin(), v.end(), v.begin(), (int(*)(int)) tolower);
                if (values.count(v) > 0) {
                    // write ATTRVALUE
                    writeByte(values[v]);
                } else {
                    if (stringValues.count(value) <= 0) {
                        // define this new string
                        writeStringDefine(value);
                        int index = stringValues.size();
                        stringValues[value] = index;
                    } else {
                        // write string ref
                        writeByte(STR_R);
                        writeMultiByteInt(stringValues[value]);    
                    }
                }
            }
    	    stringAttributes.clear();   	     

            // write color attributes
	        for (std::map<std::string,std::vector<double> >::iterator i = colorAttributes.begin(); i != colorAttributes.end(); i++) {
    		    std::string name = i->first;
    		    std::vector<double> value = i->second;
                // write ATTRSTART
                writeByte(attributes[name]);
                // write OPAQUE
                writeByte(OPAQUE);
                writeMultiByteInt(value.size());
                writeByte((int)(value[0] * 0xff) & 0xff);
                writeByte((int)(value[1] * 0xff) & 0xff);
                writeByte((int)(value[2] * 0xff) & 0xff);
                if (value.size() > 3) writeByte((int)(value[3] * 0xff) & 0xff);
            }
    	    colorAttributes.clear();
    	    
            // write long attributes
	        for (std::map<std::string,int64>::iterator i = longAttributes.begin(); i != longAttributes.end(); i++) {
    		    std::string name = i->first;
    		    int64 value = i->second;
                 // write ATTRSTART
                writeByte(attributes[name]);
                // write OPAQUE
                writeByte(OPAQUE);
                writeMultiByteInt(8);
                writeLong(value);
            }
    	    longAttributes.clear();
    	    
            // write int attributes
	        for (std::map<std::string,int>::iterator i = intAttributes.begin(); i != intAttributes.end(); i++) {
    		    std::string name = i->first;
    		    int value = i->second;
                // write ATTRSTART
                writeByte(attributes[name]);
                // write OPAQUE
                writeByte(OPAQUE);
                writeMultiByteInt(4);
                writeInt(value);
            }
    	    intAttributes.clear();
    	    
            // write boolean attributes
	        for (std::map<std::string,bool>::iterator i = booleanAttributes.begin(); i != booleanAttributes.end(); i++) {
    		    std::string name = i->first;
    		    bool value = i->second;
                // write ATTRSTART
                writeByte(attributes[name]);
                // write ATTRVALUE
                writeByte(value ? values["true"] : values["false"]);
            }
    	    booleanAttributes.clear();
    	    
            // write double attributes
	        for (std::map<std::string,double>::iterator i = doubleAttributes.begin(); i != doubleAttributes.end(); i++) {
    		    std::string name = i->first;
    		    double value = i->second;
    		    if (!isPoint && (name != "x") && (name != "y") && (name != "z")) {
                    // write ATTRSTART
                    writeByte(attributes[name]);
                    // write OPAQUE
                    writeByte(OPAQUE);
                    writeMultiByteInt(singlePrecision ? 4 : 8);
                    writeReal(value);
                }        
            }
    	    doubleAttributes.clear();
    	    
    	    // end of attributes
    	    writeByte(END);   	     
	    }
        
        if (s == "point") {
            writeByte(OPAQUE);
            writeMultiByteInt(singlePrecision ?  12 : 24);
            writeReal(doubleAttributes["x"]);
            writeReal(doubleAttributes["y"]);
            writeReal(doubleAttributes["z"]);
        }
        
        if (isPoint && !hasContent) {
            // end this tag
            writeByte(END);
        }    
    }
    
    void BHepRepWriter::writePoints() {
        if (points.size() <= 0) return;
        
        writeByte(tags["point"] | CONTENT);                
        writeByte(OPAQUE);
        writeMultiByteInt(points.size()*(singlePrecision ? 4 : 8));
        for (std::vector<double>::iterator i = points.begin(); i != points.end(); ) {
            writeReal(*i++);
            writeReal(*i++);
            writeReal(*i++);
        }
        writeByte(END);
        
        points.clear();
    }
        
    void BHepRepWriter::setAttribute(std::string name, char* value) {
        setAttribute(name, (std::string)value);
    }

    void BHepRepWriter::setAttribute(std::string name, std::string value) {
        if (name == "value") name = name.append("String");
                
        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        stringAttributes[name] = value;
    }

    void BHepRepWriter::setAttribute(std::string name, std::vector<double> value) {
        if (name == "value") name = name.append("Color");

        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        colorAttributes[name] = value;
    }
    
    void BHepRepWriter::setAttribute(std::string name, int64 value) {
        if (name == "value") name = name.append("Long");

        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        longAttributes[name] = value;
    }
    
    void BHepRepWriter::setAttribute(std::string name, int value) {
        if (name == "value") name = name.append("Int");

        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        intAttributes[name] = value;
    }
    
    void BHepRepWriter::setAttribute(std::string name, bool value) {
        if (name == "value") name = name.append("Boolean");

        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        booleanAttributes[name] = value;
    }
    
    void BHepRepWriter::setAttribute(std::string name, double value) {
        if (name == "value") name = name.append("Double");

        // make sure the attribute name is defined
        if (attributes.count(name) <= 0) {
            std::cerr << "Cannot find attribute name '" << name << "' in attributes table, skipped." << std::endl;
            return;
        }
                    
        doubleAttributes[name] = value;
    }

    void BHepRepWriter::writeStringDefine(std::string s) {
        writeByte(STR_D);
        writeString(s);
    }

    void BHepRepWriter::writeMultiByteInt(unsigned int ui) {
        unsigned char buf[5];
        int idx = 0;
        
        do {
            buf[idx++] = (unsigned char) (ui & 0x7f);
            ui = ui >> 7;
        }
        while (ui != 0);

        while (idx > 1) {
            writeByte(buf[--idx] | 0x80);
        }
        writeByte(buf[0]);
    }

    void BHepRepWriter::writeReal(double d) {
        if (singlePrecision) {
            union {
    	        int i;
    	        float f;
            } u;
            u.f = (float)d;
            writeInt(u.i);
        } else {
            union {
    	        int64 i;
    	        double d;
            } u;
            u.d = d;
     
            writeLong(u.i);
        }
    }

    void BHepRepWriter::writeLong(int64 i) {        
        // write network-order
        os.put((i >> 56) & 0xff);
        os.put((i >> 48) & 0xff);
        os.put((i >> 40) & 0xff);
        os.put((i >> 32) & 0xff);
        os.put((i >> 24) & 0xff);
        os.put((i >> 16) & 0xff);
        os.put((i >>  8) & 0xff);
        os.put((i      ) & 0xff);
    }    
    
    void BHepRepWriter::writeInt(int i) {
        // write network-order
        os.put((i >> 24) & 0xff);
        os.put((i >> 16) & 0xff);
        os.put((i >>  8) & 0xff);
        os.put((i      ) & 0xff);
    }
    
    void BHepRepWriter::writeByte(unsigned char b) {
        os.put((char)b);
    }
    
    void BHepRepWriter::writeString(std::string s) {
        os << s;
        os.put(0);
    }
} // cheprep
