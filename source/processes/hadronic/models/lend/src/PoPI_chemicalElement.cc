/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <ctype.h>

#include <PoPI.hpp>

namespace PoPI {

#define PoPI_chemicalElementChars "chemicalElement"
#define PoPI_isotopesChars "isotopes"
#define PoPI_Z_Chars "Z"

static std::string emptyString( "" );
static std::map<int, std::string>  ZtoChemicalElementSymbols{
    {1,   "H"},  {2,   "He"}, {3,   "Li"}, {4,   "Be"}, {5,   "B"},  {6,   "C"},  {7,   "N"},  {8,   "O"},  {9,   "F"},  {10,  "Ne"},
    {11,  "Na"}, {12,  "Mg"}, {13,  "Al"}, {14,  "Si"}, {15,  "P"},  {16,  "S"},  {17,  "Cl"}, {18,  "Ar"}, {19,  "K"},  {20,  "Ca"},
    {21,  "Sc"}, {22,  "Ti"}, {23,  "V"},  {24,  "Cr"}, {25,  "Mn"}, {26,  "Fe"}, {27,  "Co"}, {28,  "Ni"}, {29,  "Cu"}, {30,  "Zn"},
    {31,  "Ga"}, {32,  "Ge"}, {33,  "As"}, {34,  "Se"}, {35,  "Br"}, {36,  "Kr"}, {37,  "Rb"}, {38,  "Sr"}, {39,  "Y"},  {40,  "Zr"},
    {41,  "Nb"}, {42,  "Mo"}, {43,  "Tc"}, {44,  "Ru"}, {45,  "Rh"}, {46,  "Pd"}, {47,  "Ag"}, {48,  "Cd"}, {49,  "In"}, {50,  "Sn"},
    {51,  "Sb"}, {52,  "Te"}, {53,  "I"},  {54,  "Xe"}, {55,  "Cs"}, {56,  "Ba"}, {57,  "La"}, {58,  "Ce"}, {59,  "Pr"}, {60,  "Nd"},
    {61,  "Pm"}, {62,  "Sm"}, {63,  "Eu"}, {64,  "Gd"}, {65,  "Tb"}, {66,  "Dy"}, {67,  "Ho"}, {68,  "Er"}, {69,  "Tm"}, {70,  "Yb"},
    {71,  "Lu"}, {72,  "Hf"}, {73,  "Ta"}, {74,  "W"},  {75,  "Re"}, {76,  "Os"}, {77,  "Ir"}, {78,  "Pt"}, {79,  "Au"}, {80,  "Hg"},
    {81,  "Tl"}, {82,  "Pb"}, {83,  "Bi"}, {84,  "Po"}, {85,  "At"}, {86,  "Rn"}, {87,  "Fr"}, {88,  "Ra"}, {89,  "Ac"}, {90,  "Th"},
    {91,  "Pa"}, {92,  "U"},  {93,  "Np"}, {94,  "Pu"}, {95,  "Am"}, {96,  "Cm"}, {97,  "Bk"}, {98,  "Cf"}, {99,  "Es"}, {100, "Fm"},
    {101, "Md"}, {102, "No"}, {103, "Lr"}, {104, "Rf"}, {105, "Db"}, {106, "Sg"}, {107, "Bh"}, {108, "Hs"}, {109, "Mt"}, {110, "Ds"},
    {111, "Rg"}, {112, "Cn"}, {113, "Nh"}, {114, "Fl"}, {115, "Mc"}, {116, "Lv"}, {117, "Ts"}, {118, "Og"} };

static std::map<std::string, int> chemicalElementSymbolToZs;
std::map<std::string, std::string> supportedNucleusAliases{ {"d", "h2"}, {"t", "h3"}, {"h", "he3"}, {"a", "he4"} };
static std::string protonFakeAlias( "h1" );

/*! \class ChemicalElement
 * This class represents a **PoPs** chemicalElement instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **GNDS** chemicalElement node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **ChemicalElement** to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

ChemicalElement::ChemicalElement( HAPI::Node const &a_node, Database *a_DB, LUPI_maybeUnused Database *a_parent ) :
        SymbolBase( a_node, Particle_class::chemicalElement ),
        m_Z( a_node.attribute( PoPI_Z_Chars ).as_int( ) ),
        m_name( a_node.attribute( PoPI_nameChars ).value( ) ),
        m_isotopes( PoPI_isotopesChars ) {

    addToSymbols( a_DB );
    m_isotopes.appendFromParentNode( a_node.child( PoPI_isotopesChars ), a_DB, this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ChemicalElement::~ChemicalElement( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void ChemicalElement::calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const {

    for( std::size_t i1 = 0; i1 <  m_isotopes.size( ); ++i1 ) {
        Isotope const &isotope = m_isotopes[i1];

        isotope.calculateNuclideGammaBranchStateInfos( a_pops, a_nuclideGammaBranchStateInfos );
    }
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void ChemicalElement::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string::size_type size = m_isotopes.size( );
    std::string ZStr = LUPI::Misc::argumentsToString( "%d", m_Z );

    if( size == 0 ) return;

    std::string header = a_indent1 + "<chemicalElement symbol=\"" + symbol( ) + "\" Z=\"" + ZStr + "\" name=\"" + m_name + "\">";
    a_XMLList.push_back( header );

    std::string indent2 = a_indent1 + "  ";
    std::string isotopeSuite = indent2 + "<" + PoPI_isotopesChars + ">";
    a_XMLList.push_back( isotopeSuite );

    std::string indent3 = indent2 + "  ";
    for( std::string::size_type i1 = 0; i1 < size; ++i1 ) m_isotopes[i1].toXMLList( a_XMLList, indent3 );

    appendXMLEnd( a_XMLList, PoPI_isotopesChars );
    appendXMLEnd( a_XMLList, PoPI_chemicalElementChars );
}

/* *********************************************************************************************************//**
 * Returns the maximum supported Z (atomic number) supported by function chemicalElementInfoFromZ.
 *
 * @return                          int.
 ***********************************************************************************************************/

int maximumChemicalElementZ( ) {

    return( static_cast<int>( ZtoChemicalElementSymbols.size( ) ) );
}

/* *********************************************************************************************************//**
 * Returns the chemical element's symbol (*a_wantSymbol* = **true**) or name (*a_wantSymbol* = **false**) for
 * the requrest atomic number *a_Z*. *a_wantSymbol* is **true**, the returned symbol is for the nuclide when
 * *a_asNucleus* is **false** and for the nucleus otherwise.
 *
 * @param a_Z               [in]    The Z (atomic number of the chemical element.
 * @param a_wantSymbol      [in]    If **true** returns the chemical element's symbol otherwise its name.
 * @param a_asNucleus       [in]    If **true** returns the symbol for the nucleus, otherwise for the nuclide. As no affect if a_wantSymbol is **false**.
 *
 * @return                          The symbol for the nuclide or the nucleus, or the name as a std::string.
 ***********************************************************************************************************/

std::string chemicalElementInfoFromZ( int a_Z, bool a_wantSymbol, bool a_asNucleus ) {

    std::string info;

    switch( a_Z ) {
    case 1:
        if( a_wantSymbol ) {
            info = "H"; }
        else {
            info = "Hydrogen";
        }
        break;
    case 2:
        if( a_wantSymbol ) {
            info = "He"; }
        else {
            info = "Helium";
        }
        break;
    case 3:
        if( a_wantSymbol ) {
            info = "Li"; }
        else {
            info = "Lithium";
        }
        break;
    case 4:
        if( a_wantSymbol ) {
            info = "Be"; }
        else {
            info = "Beryllium";
        }
        break;
    case 5:
        if( a_wantSymbol ) {
            info = "B"; }
        else {
            info = "Boron";
        }
        break;
    case 6:
        if( a_wantSymbol ) {
            info = "C"; }
        else {
            info = "Carbon";
        }
        break;
    case 7:
        if( a_wantSymbol ) {
            info = "N"; }
        else {
            info = "Nitrogen";
        }
        break;
    case 8:
        if( a_wantSymbol ) {
            info = "O"; }
        else {
            info = "Oxygen";
        }
        break;
    case 9:
        if( a_wantSymbol ) {
            info = "F"; }
        else {
            info = "Fluorine";
        }
        break;
    case 10:
        if( a_wantSymbol ) {
            info = "Ne"; }
        else {
            info = "Neon";
        }
        break;
    case 11:
        if( a_wantSymbol ) {
            info = "Na"; }
        else {
            info = "Sodium";
        }
        break;
    case 12:
        if( a_wantSymbol ) {
            info = "Mg"; }
        else {
            info = "Magnesium";
        }
        break;
    case 13:
        if( a_wantSymbol ) {
            info = "Al"; }
        else {
            info = "Aluminium";
        }
        break;
    case 14:
        if( a_wantSymbol ) {
            info = "Si"; }
        else {
            info = "Silicon";
        }
        break;
    case 15:
        if( a_wantSymbol ) {
            info = "P"; }
        else {
            info = "Phosphorus";
        }
        break;
    case 16:
        if( a_wantSymbol ) {
            info = "S"; }
        else {
            info = "Sulphur";
        }
        break;
    case 17:
        if( a_wantSymbol ) {
            info = "Cl"; }
        else {
            info = "Chlorine";
        }
        break;
    case 18:
        if( a_wantSymbol ) {
            info = "Ar"; }
        else {
            info = "Argon";
        }
        break;
    case 19:
        if( a_wantSymbol ) {
            info = "K"; }
        else {
            info = "Potassium";
        }
        break;
    case 20:
        if( a_wantSymbol ) {
            info = "Ca"; }
        else {
            info = "Calcium";
        }
        break;
    case 21:
        if( a_wantSymbol ) {
            info = "Sc"; }
        else {
            info = "Scandium";
        }
        break;
    case 22:
        if( a_wantSymbol ) {
            info = "Ti"; }
        else {
            info = "Titanium";
        }
        break;
    case 23:
        if( a_wantSymbol ) {
            info = "V"; }
        else {
            info = "Vanadium";
        }
        break;
    case 24:
        if( a_wantSymbol ) {
            info = "Cr"; }
        else {
            info = "Chromium";
        }
        break;
    case 25:
        if( a_wantSymbol ) {
            info = "Mn"; }
        else {
            info = "Manganese";
        }
        break;
    case 26:
        if( a_wantSymbol ) {
            info = "Fe"; }
        else {
            info = "Iron";
        }
        break;
    case 27:
        if( a_wantSymbol ) {
            info = "Co"; }
        else {
            info = "Cobalt";
        }
        break;
    case 28:
        if( a_wantSymbol ) {
            info = "Ni"; }
        else {
            info = "Nickel";
        }
        break;
    case 29:
        if( a_wantSymbol ) {
            info = "Cu"; }
        else {
            info = "Copper";
        }
        break;
    case 30:
        if( a_wantSymbol ) {
            info = "Zn"; }
        else {
            info = "Zinc";
        }
        break;
    case 31:

        if( a_wantSymbol ) {
            info = "Ga"; }
        else {
            info = "Gallium";
        }
        break;
    case 32:
        if( a_wantSymbol ) {
            info = "Ge"; }
        else {
            info = "Germanium";
        }
        break;
    case 33:
        if( a_wantSymbol ) {
            info = "As"; }
        else {
            info = "Arsenic";
        }
        break;
    case 34:
        if( a_wantSymbol ) {
            info = "Se"; }
        else {
            info = "Selenium";
        }
        break;
    case 35:
        if( a_wantSymbol ) {
            info = "Br"; }
        else {
            info = "Bromine";
        }
        break;
    case 36:
        if( a_wantSymbol ) {
            info = "Kr"; }
        else {
            info = "Krypton";
        }
        break;
    case 37:
        if( a_wantSymbol ) {
            info = "Rb"; }
        else {
            info = "Rubidium";
        }
        break;
    case 38:
        if( a_wantSymbol ) {
            info = "Sr"; }
        else {
            info = "Strontium";
        }
        break;
    case 39:
        if( a_wantSymbol ) {
            info = "Y"; }
        else {
            info = "Yttrium";
        }
        break;
    case 40:
        if( a_wantSymbol ) {
            info = "Zr"; }
        else {
            info = "Zirconium";
        }
        break;
    case 41:
        if( a_wantSymbol ) {
            info = "Nb"; }
        else {
            info = "Niobium";
        }
        break;
    case 42:
        if( a_wantSymbol ) {
            info = "Mo"; }
        else {
            info = "Molybdenum";
        }
        break;
    case 43:
        if( a_wantSymbol ) {
            info = "Tc"; }
        else {
            info = "Technetium";
        }
        break;
    case 44:
        if( a_wantSymbol ) {
            info = "Ru"; }
        else {
            info = "Ruthenium";
        }
        break;
    case 45:
        if( a_wantSymbol ) {
            info = "Rh"; }
        else {
            info = "Rhodium";
        }
        break;
    case 46:
        if( a_wantSymbol ) {
            info = "Pd"; }
        else {
            info = "Palladium";
        }
        break;
    case 47:
        if( a_wantSymbol ) {
            info = "Ag"; }
        else {
            info = "Silver";
        }
        break;
    case 48:
        if( a_wantSymbol ) {
            info = "Cd"; }
        else {
            info = "Cadmium";
        }
        break;
    case 49:
        if( a_wantSymbol ) {
            info = "In"; }
        else {
            info = "Indium";
        }
        break;
    case 50:
        if( a_wantSymbol ) {
            info = "Sn"; }
        else {
            info = "Tin";
        }
        break;
    case 51:
        if( a_wantSymbol ) {
            info = "Sb"; }
        else {
            info = "Antimony";
        }
        break;
    case 52:
        if( a_wantSymbol ) {
            info = "Te"; }
        else {
            info = "Tellurium";
        }
        break;
    case 53:
        if( a_wantSymbol ) {
            info = "I"; }
        else {
            info = "Iodine";
        }
        break;
    case 54:
        if( a_wantSymbol ) {
            info = "Xe"; }
        else {
            info = "Xenon";
        }
        break;
    case 55:
        if( a_wantSymbol ) {
            info = "Cs"; }
        else {
            info = "Cesium";
        }
        break;
    case 56:
        if( a_wantSymbol ) {
            info = "Ba"; }
        else {
            info = "Barium";
        }
        break;
    case 57:
        if( a_wantSymbol ) {
            info = "La"; }
        else {
            info = "Lanthanum";
        }
        break;
    case 58:
        if( a_wantSymbol ) {
            info = "Ce"; }
        else {
            info = "Cerium";
        }
        break;
    case 59:
        if( a_wantSymbol ) {
            info = "Pr"; }
        else {
            info = "Praseodymium";
        }
        break;
    case 60:
        if( a_wantSymbol ) {
            info = "Nd"; }
        else {
            info = "Neodymium";
        }
        break;
    case 61:
        if( a_wantSymbol ) {
            info = "Pm"; }
        else {
            info = "Promethium";
        }
        break;
    case 62:
        if( a_wantSymbol ) {
            info = "Sm"; }
        else {
            info = "Samarium";
        }
        break;
    case 63:
        if( a_wantSymbol ) {
            info = "Eu"; }
        else {
            info = "Europium";
        }
        break;
    case 64:
        if( a_wantSymbol ) {
            info = "Gd"; }
        else {
            info = "Gadolinium";
        }
        break;
    case 65:
        if( a_wantSymbol ) {
            info = "Tb"; }
        else {
            info = "Terbium";
        }
        break;
    case 66:
        if( a_wantSymbol ) {
            info = "Dy"; }
        else {
            info = "Dysprosium";
        }
        break;
    case 67:
        if( a_wantSymbol ) {
            info = "Ho"; }
        else {
            info = "Holmium";
        }
        break;
    case 68:
        if( a_wantSymbol ) {
            info = "Er"; }
        else {
            info = "Erbium";
        }
        break;
    case 69:
        if( a_wantSymbol ) {
            info = "Tm"; }
        else {
            info = "Thulium";
        }
        break;
    case 70:
        if( a_wantSymbol ) {
            info = "Yb"; }
        else {
            info = "Ytterbium";
        }
        break;
    case 71:
        if( a_wantSymbol ) {
            info = "Lu"; }
        else {
            info = "Lutetium";
        }
        break;
    case 72:
        if( a_wantSymbol ) {
            info = "Hf"; }
        else {
            info = "Hafnium";
        }
        break;
    case 73:
        if( a_wantSymbol ) {
            info = "Ta"; }
        else {
            info = "Tantalum";
        }
        break;
    case 74:
        if( a_wantSymbol ) {
            info = "W"; }
        else {
            info = "Tungsten";
        }
        break;
    case 75:
        if( a_wantSymbol ) {
            info = "Re"; }
        else {
            info = "Rhenium";
        }
        break;
    case 76:
        if( a_wantSymbol ) {
            info = "Os"; }
        else {
            info = "Osmium";
        }
        break;
    case 77:
        if( a_wantSymbol ) {
            info = "Ir"; }
        else {
            info = "Iridium";
        }
        break;
    case 78:
        if( a_wantSymbol ) {
            info = "Pt"; }
        else {
            info = "Platinum";
        }
        break;
    case 79:
        if( a_wantSymbol ) {
            info = "Au"; }
        else {
            info = "Gold";
        }
        break;
    case 80:
        if( a_wantSymbol ) {
            info = "Hg"; }
        else {
            info = "Mercury";
        }
        break;
    case 81:
        if( a_wantSymbol ) {
            info = "Tl"; }
        else {
            info = "Thallium";
        }
        break;
    case 82:
        if( a_wantSymbol ) {
            info = "Pb"; }
        else {
            info = "Lead";
        }
        break;
    case 83:
        if( a_wantSymbol ) {
            info = "Bi"; }
        else {
            info = "Bismuth";
        }
        break;
    case 84:
        if( a_wantSymbol ) {
            info = "Po"; }
        else {
            info = "Polonium";
        }
        break;
    case 85:
        if( a_wantSymbol ) {
            info = "At"; }
        else {
            info = "Astatine";
        }
        break;
    case 86:
        if( a_wantSymbol ) {
            info = "Rn"; }
        else {
            info = "Radon";
        }
        break;
    case 87:
        if( a_wantSymbol ) {
            info = "Fr"; }
        else {
            info = "Francium";
        }
        break;
    case 88:
        if( a_wantSymbol ) {
            info = "Ra"; }
        else {
            info = "Radium";
        }
        break;
    case 89:
        if( a_wantSymbol ) {
            info = "Ac"; }
        else {
            info = "Actinium";
        }
        break;
    case 90:
        if( a_wantSymbol ) {
            info = "Th"; }
        else {
            info = "Thorium";
        }
        break;
    case 91:
        if( a_wantSymbol ) {
            info = "Pa"; }
        else {
            info = "Protactinium";
        }
        break;
    case 92:
        if( a_wantSymbol ) {
            info = "U"; }
        else {
            info = "Uranium";
        }
        break;
    case 93:
        if( a_wantSymbol ) {
            info = "Np"; }
        else {
            info = "Neptunium";
        }
        break;
    case 94:
        if( a_wantSymbol ) {
            info = "Pu"; }
        else {
            info = "Plutonium";
        }
        break;
    case 95:
        if( a_wantSymbol ) {
            info = "Am"; }
        else {
            info = "Americium";
        }
        break;
    case 96:
        if( a_wantSymbol ) {
            info = "Cm"; }
        else {
            info = "Curium";
        }
        break;
    case 97:
        if( a_wantSymbol ) {
            info = "Bk"; }
        else {
            info = "Berkelium";
        }
        break;
    case 98:
        if( a_wantSymbol ) {
            info = "Cf"; }
        else {
            info = "Californium";
        }
        break;
    case 99:
        if( a_wantSymbol ) {
            info = "Es"; }
        else {
            info = "Einsteinium";
        }
        break;
    case 100:
        if( a_wantSymbol ) {
            info = "Fm"; }
        else {
            info = "Fermium";
        }
        break;
    case 101:
        if( a_wantSymbol ) {
            info = "Md"; }
        else {
            info = "Mendelevium";
        }
        break;
    case 102:
        if( a_wantSymbol ) {
            info = "No"; }
        else {
            info = "Nobelium";
        }
        break;
    case 103:
        if( a_wantSymbol ) {
            info = "Lr"; }
        else {
            info = "Lawrencium";
        }
        break;
    case 104:
        if( a_wantSymbol ) {
            info = "Rf"; }
        else {
            info = "Rutherfordium";
        }
        break;
    case 105:
        if( a_wantSymbol ) {
            info = "Db"; }
        else {
            info = "Dubnium";
        }
        break;
    case 106:
        if( a_wantSymbol ) {
            info = "Sg"; }
        else {
            info = "Seaborgium";
        }
        break;
    case 107:
        if( a_wantSymbol ) {
            info = "Bh"; }
        else {
            info = "Bohrium";
        }
        break;
    case 108:
        if( a_wantSymbol ) {
            info = "Hs"; }
        else {
            info = "Hassium";
        }
        break;
    case 109:
        if( a_wantSymbol ) {
            info = "Mt"; }
        else {
            info = "Meitnerium";
        }
        break;
    case 110:
        if( a_wantSymbol ) {
            info = "Ds"; }
        else {
            info = "Darmstadtium";
        }
        break;
    case 111:
        if( a_wantSymbol ) {
            info = "Rg"; }
        else {
            info = "Roentgenium";
        }
        break;
    case 112:
        if( a_wantSymbol ) {
            info = "Cn"; }
        else {
            info = "Copernicium";
        }
        break;
    case 113:
        if( a_wantSymbol ) {
            info = "Nh"; }
        else {
            info = "Nihonium";
        }
        break;
    case 114:
        if( a_wantSymbol ) {
            info = "Fl"; }
        else {
            info = "Flerovium";
        }
        break;
    case 115:
        if( a_wantSymbol ) {
            info = "Mc"; }
        else {
            info = "Moscovium";
        }
        break;
    case 116:
        if( a_wantSymbol ) {
            info = "Lv"; }
        else {
            info = "Livermorium";
        }
        break;
    case 117:
        if( a_wantSymbol ) {
            info = "Ts"; }
        else {
            info = "Tennessine";
        }
        break;
    case 118:
        if( a_wantSymbol ) {
            info = "Og"; }
        else {
            info = "Oganesson";
        }
        break;
    default:
        break;
    }

    if( a_wantSymbol && a_asNucleus ) {
        char c1[3];
        c1[0] = tolower( info.c_str( )[0] );
        c1[1] = 0;
        c1[2] = 0;
        if( info.size( ) > 1 ) c1[1] = info.c_str( )[1];

        info = c1;
    }

    return( info );
}

/* *********************************************************************************************************//**
 * Returns the chemical element symbol for the requested atomic number *a_Z*.
 *
 * @param a_Z               [in]    The atomic number (Z) of the requested chemical element.
 *
 * @return                          The symbol for the nuclide or an empty string if *a_Z* is an invalid atomic number.
 ***********************************************************************************************************/

std::string const &chemicalElementSymbolFromZ( int a_Z ) {

    if( ZtoChemicalElementSymbols.find( a_Z ) == ZtoChemicalElementSymbols.end( ) ) return( emptyString );

    return( ZtoChemicalElementSymbols[a_Z ] );
}

/* *********************************************************************************************************//**
 * Returns the atomic number (Z) for the requested chemical element's symbol.
 *
 * @param a_symbol          [in]    The atomic symbol.
 *
 * @return                          The atomic number or 0 if *a_symbol* is an invalid symbol.
 ***********************************************************************************************************/

int Z_FromChemicalElementSymbol( std::string const &a_symbol ) {

    if( chemicalElementSymbolToZs.size( ) == 0 ) {
        for( auto iter = ZtoChemicalElementSymbols.begin( ); iter != ZtoChemicalElementSymbols.end( ); ++iter ) {
            chemicalElementSymbolToZs[iter->second] = iter->first;
        }
    }

    if( chemicalElementSymbolToZs.find( a_symbol ) == chemicalElementSymbolToZs.end( ) ) return( 0 );

    return( chemicalElementSymbolToZs[a_symbol] );
}


/* *********************************************************************************************************//**
 * This class breaks down a PoPs id for a nuclide or nuclear into its components (e.g., Z, A, index).
 * The *a_id* can also be a nuclear meta-stable or one of the light paritlce aliases (i.e., "d", "t", "h" or "a").
 * If *a_id* is a light particle alias, its nucleus equavalent is used. Also, "p" is treated as "h1", and "n"
 * returns Z = 0 and A = 1. Currently, no other PoPs id's are supported.
 *
 * @param a_id              [in]    The PoPs id of the particle.
 ***********************************************************************************************************/

ParseIdInfo::ParseIdInfo( std::string const &a_id ) :
        m_isSupported( false ),
        m_id( a_id ),
        m_isNuclear( false ),
        m_isNucleus( false ),
        m_isChemicalElement( false ),
        m_isAnti( false ),
        m_isMetaStable( false ),
        m_symbol( "" ),
        m_Z( 0 ),
        m_A( 0 ), 
        m_index( 0 ),
        m_qualifier( "" ) {

    std::string a_anti;

    std::string baseId = baseAntiQualifierFromID( a_id, a_anti, &m_qualifier );
    m_isAnti = IDs::anti == a_anti;

    if( supportedNucleusAliases.find( baseId ) != supportedNucleusAliases.end( ) ) {
        baseId = supportedNucleusAliases[baseId]; }
    else if( baseId == IDs::proton ) {
        baseId = protonFakeAlias;
    }

    if( baseId == "n" ) {
        m_A = 1;
        m_isSupported = true;
        return;
    }

    std::vector<std::string> parts;
    if( baseId.find( "_m" ) != std::string::npos ) {
        m_isMetaStable = true;
        parts = LUPI::Misc::splitString( baseId, "_m" ); }
    else {
        parts = LUPI::Misc::splitString( baseId, "_e" );
    }

// Now look for something of the form "SA(_[em]N)" in parts[0] where S is symbol, A is atomic number and "_[em]N" is options nulcear level or
// meta-stable specifier. If no match is found, assume a_id does not define a nuclear id.
    std::string isotope = parts[0];
    std::size_t digitIndex = isotope.find_first_of( "01233456789" );

    std::string symbol( isotope.substr( 0, digitIndex ) );                      // This should be S.
    std::string symbolCap;
    if( symbol.size( ) > 0 ) {
        char firstChar[2];
        firstChar[0] = std::toupper( symbol[0] );
        firstChar[1] = 0;
        std::string firstStringChar( firstChar );
        symbolCap = firstStringChar + symbol.substr( 1 );
    }

    if( digitIndex != std::string::npos ) {
        std::string AStr( isotope.substr( digitIndex ) );                       // This should be A.
        if( symbol.size( ) > 0 ) {
            m_Z = Z_FromChemicalElementSymbol( symbolCap );
            if( m_Z > 0 ) {                                                     // We have a valid chemical element symbol.
                if( ( AStr.size( ) > 0 ) && ( LUPI::Misc::stringToInt( AStr, m_A ) ) ) {
                    if( m_A < 0 ) {
                        m_A = 0; }
                    else {
                        bool isValidNuclearId = parts.size( ) == 1;

                        if( parts.size( ) > 1 ) {
                            isValidNuclearId = ( parts.size( ) == 2 ) && LUPI::Misc::stringToInt( parts[1], m_index );
                        }

                        if( isValidNuclearId ) {                                // Should be a valid nuclear id.
                            m_symbol = symbolCap;
                            m_isNuclear = true;
                            m_isNucleus =  symbolCap != symbol;
                        }
                    }
                }
                m_isSupported = true;
            }
        } }
    else if( symbol.size( ) > 0 ) {
        m_Z = Z_FromChemicalElementSymbol( symbolCap );
        if( m_Z > 0 ) {
            m_symbol = symbolCap;
            m_isChemicalElement = true;
            m_isSupported = true;
        }
    }
}

/* *********************************************************************************************************//**
 * This method prints the contents of *this*. This is mainly for debugging.
 *
 * @param a_terse           [in]    If **true**, all members are printed on one line with no description. Otherwise, each member is printed on a separate line with a description.
 * @param a_indent          [in]    The amount of indentation on each line before anything is printed.
 ***********************************************************************************************************/

void ParseIdInfo::print( bool a_terse, std::string const &a_indent ) const {

    if( a_terse ) {
        std::cout << a_indent << m_id 
                << boolToString( m_isSupported, " " ).c_str( ) 
                << boolToString( m_isNuclear, " " ).c_str( ) 
                << boolToString( m_isNucleus, " " ).c_str( )
                << boolToString( m_isChemicalElement, " " ).c_str( ) 
                << boolToString( m_isAnti, " " ).c_str( ) 
                << boolToString( m_isMetaStable, " " ).c_str( ) 
                << LUPI::Misc::argumentsToString( " %s", m_symbol.c_str( ) ) 
                << LUPI::Misc::argumentsToString( " %d", m_Z ) 
                << LUPI::Misc::argumentsToString( " %d", m_A ) 
                << LUPI::Misc::argumentsToString( " %d", m_index )
                << LUPI::Misc::argumentsToString( " %s", m_qualifier.c_str( ) ) 
                << std::endl; }
    else {
        std::cout << a_indent << "id = " << m_id << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isSupported = %s", boolToString( m_isSupported, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isNuclear = %s", boolToString( m_isNuclear, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isNucleus = %s", boolToString( m_isNucleus, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isChemicalElement = %s", boolToString( m_isChemicalElement, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isAnti = %s", boolToString( m_isAnti, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  isMetaStable = %s", boolToString( m_isMetaStable, "" ).c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  symbol = <%s>", m_symbol.c_str( ) ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  Z = %d", m_Z ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  A = %d", m_A ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  index = %d", m_index ) << std::endl;
        std::cout << a_indent << LUPI::Misc::argumentsToString( "  qualifier = <%s>", m_qualifier.c_str( ) ) << std::endl;
    }
}

/* *********************************************************************************************************//**
 * This method returns a string representation of *a_value*.
 *
 * @param a_value           [in]    If **true**, all members are printed on one line with no description. Otherwise, each member is printed on a separate line w
 * @param a_prefix          [in]    The amount of indentation on each line before anything is printed.
 *
 * @return                          The string representation of *a_value*.
 ***********************************************************************************************************/

std::string ParseIdInfo::boolToString( bool a_value, std::string const &a_prefix ) const {

    std::string boolString( a_prefix );

    if( a_value ) {
        boolString += "true"; }
    else {
        boolString += "false";
    }

    return( boolString );
}

}
