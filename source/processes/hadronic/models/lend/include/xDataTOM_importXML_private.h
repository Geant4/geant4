/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef xDataTOM_importXML_private_h_included
#define xDataTOM_importXML_private_h_included

#include <expat.h>
#include <statusMessageReporting.h>
#include <xDataTOM_private.h>

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

#ifndef XMLCALL
#define XMLCALL
#endif

#ifndef XML_Char
#define XML_Char char
#endif

#ifndef XML_Size
#define XML_Size long
#endif

typedef struct xDataXMLType_s xDataXMLType;

enum xDataXML_status { xDataXML_statusParsing = 1, xDataXML_statusCompleted, xDataXML_statusError };
enum xDataXML_errorCodes { xDataXML_errNone, xDataXML_errXML_ParserCreate, xDataXML_errFileError, xDataXML_errXMLParser };
enum xDataXML_itemMode { xDataXML_itemModeEnd = 0, xDataXML_itemModeElement, xDataXML_itemModeText };
enum xDataXML_xDataType { xDataXML_xDataType_Ok, xDataXML_xDataType_ConvertingToData, xDataXML_xDataType_ConvertingToString };

typedef struct xDataTOM_importXMLTOM_s xDataTOM_importXMLTOM;
typedef struct xDataXML_attribute_s xDataXML_attribute;
typedef struct xDataXML_document_s xDataXML_document;
typedef struct xDataXML_docInfo_s xDataXML_docInfo;
typedef struct xDataXML_attributionList_s xDataXML_attributionList;
typedef struct xDataXML_element_s xDataXML_element;
typedef struct xDataXML_elementList_s xDataXML_elementList;
typedef struct xDataXML_elementListItem_s xDataXML_elementListItem;
typedef struct xDataXML_rootElement_s xDataXML_rootElement;
typedef struct xDataXML_text_s xDataXML_text;
typedef struct xDataXML_smr_s xDataXML_smr;
typedef struct xDataXML_item_s xDataXML_item;

struct xDataXML_docInfo_s {
    size_t line;
    size_t column;
};

struct xDataXML_attribute_s {
    char *name;
    char *value;
};

struct xDataXML_attributionList_s {
    int number;
    size_t size;
    xDataXML_attribute *attributes;
};

struct xDataXML_text_s {
    xDataXML_docInfo docInfo;
    size_t allocated;
    size_t length;
    char *text;
};

typedef int (*xDTXML_toDataFunction)( statusMessageReporting *smr, xDataXMLType *, xDataXML_attributionList *, char const * );
typedef char *(*xDTXML_toStringFunction)( statusMessageReporting *smr, xDataXMLType * );
typedef int (*xDTXML_releaseFunction)( statusMessageReporting *smr, xDataXMLType * );

struct xDataXML_elementListItem_s {
    xDataXML_element *element;
    const char *sortString;
};

struct xDataXML_elementList_s {
    int n;
    xDataXML_elementListItem *items;
};

struct xDataXMLType_s {
    enum xDataXML_xDataType status;
    const char *ID;
    xDataXML_element *element;
    xDTXML_toDataFunction toData;
    xDTXML_toStringFunction toString;
    xDTXML_releaseFunction release;
    int indexPresent, startPresent, endPresent, lengthPresent;
    xDataTOM_Int index, start, end, length;
    void *data;
};

struct xDataXML_rootElement_s {
    xDataXML_document *xData_doc;
    xDataXML_element *parentElement;
    xDataXML_rootElement *parentRoot;
    int depth;
    int numberOfElements;
    xDataXML_element *children;
    xDataXML_element *currentChild;
};

struct xDataXML_element_s {
    xDataXML_docInfo docInfo;
    int ordinal;                                    /* Counting from 0. */
    int index;                                      /* Value from "index" attribute if present or -1 */
    int accessed;                                   /* For the convenience of the users, not used internally. */
    xDataXML_rootElement *parentRoot;
    xDataXML_rootElement childrenRoot;
    xDataXML_element *next;
    char *name;                                     /* Allocated in xData_parseAddElementToRoot. */
    char *fullName;                                 /* Allocated in xData_parseAddElementToRoot. */
    xDataXML_attributionList attributes;               /* attributes->abbributes is allocated in xData_parseAddElementToRoot. */
    xDataXMLType xDataTypeInfo;
    size_t textOffset;
    xDataXML_text text;
};

struct xDataXML_smr_s {
    smr_userInterface smrUserInterface;
    xDataXML_document *doc;
};

struct xDataXML_document_s {
    enum xDataXML_status status;
    enum xDataXML_errorCodes error;
    enum XML_Error err;
    XML_Size err_line, err_column;
    char *fileName;
    char *realFileName;
    xDataXML_smr smrUserInterface;
    statusMessageReporting *smr;
    XML_Parser xmlParser;
    xDataXML_rootElement root;
    xDataXML_rootElement *currentRoot;
};

struct xDataXML_item_s {
    xDataXML_element *parentElement;
    xDataXML_element *element;
    enum xDataXML_itemMode mode;
    size_t textOffset;
    size_t textLength;
    char *text;
};

xDataTOM_TOM *xDataXML_importFile( statusMessageReporting *smr, char const *fileName );
xDataXML_document *xDataXML_importFile2( statusMessageReporting *smr, char const *fileName );
void *xDataXML_freeDoc( statusMessageReporting *smr, xDataXML_document *doc );
int xDataXML_parseIsError( xDataXML_document *doc );
xDataXML_element *xDataXML_getDocumentsElement( xDataXML_document *doc );
xDataXML_element *xDataXML_getFirstElement( xDataXML_element *element );
xDataXML_element *xDataXML_getNextElement( xDataXML_element *element );
enum xDataXML_itemMode xDataXML_getFirstItem( xDataXML_element *element, xDataXML_item *item );
enum xDataXML_itemMode xDataXML_getNextItem( xDataXML_item *item );
int xDataXML_isAttributeInList( xDataXML_attributionList *attributes, char const *name );
int xDataXML_isAttributeInElement( xDataXML_element *element, char const *name );
char *xDataXML_getAttributesValue( xDataXML_attributionList *attributes, char const *name );
char const *xDataXML_getAttributesValueInElement( xDataXML_element *element, char const *name );
int xDataXML_attributeListLength( xDataXML_attributionList *attributes );
xDataXML_attribute *xDataXML_attributeByIndex( xDataXML_attributionList *attributes, int index );
int xDataXML_getCommonData( statusMessageReporting *smr, xDataXML_element *element, xDataTOM_Int *index, xDataTOM_Int *start, xDataTOM_Int *end,
        xDataTOM_Int *length );
int xDataXML_xDataTypeConvertAttributes( statusMessageReporting *smr, xDataXML_element *element );
xDataTOM_Int xDataXML_convertAttributeTo_xDataTOM_Int( statusMessageReporting *smr, xDataXML_element *element, char const *name, xDataTOM_Int *n, int required );
int xDataXML_convertAttributeToDouble( statusMessageReporting *smr, xDataXML_element *element, char const *name, double *d, int required );
int xDataXML_numberOfElementsByTagName( statusMessageReporting *smr, xDataXML_element *element, char const *tagName );
xDataXML_elementList *xDataXML_getElementsByTagName( statusMessageReporting *smr, xDataXML_element *element, char const *tagName );
xDataXML_element *xDataXML_getOneElementByTagName( statusMessageReporting *smr, xDataXML_element *element, char *name, int required );
void xDataXML_freeElementList( statusMessageReporting *smr, xDataXML_elementList *list );
int xDataXML_is_xDataType( statusMessageReporting *smr, xDataXMLType *xDT, char const * const type, int setMsg );
char const *xDataXML_getFileName( xDataXML_document *doc );
char const *xDataXML_getRealFileName( xDataXML_document *doc );
xDataXML_document *xDataXML_getElementsDocument( xDataXML_element *element );
void *xDataXML_get_smrUserInterfaceFromDocument( xDataXML_document *doc );
void *xDataXML_get_smrUserInterfaceFromElement( xDataXML_element *element );
int xDataXML_stringTo_xDataTOM_Int( statusMessageReporting *smr, void *smrUserInterface, char const *c, xDataTOM_Int *value, char const *endings, char **e );
int xDataXML_stringTo_double( statusMessageReporting *smr, void *smrUserInterface, char const *c, double *value, char const *endings, char **e );
int xDataXML_addToAccessed( statusMessageReporting *smr, xDataXML_element *element, int increment );
int xDataXML_getAccessed( statusMessageReporting *smr, xDataXML_element *element );
void *xDataXML_initializeData( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE, char const *ID, size_t size );

/*
c Stuff in xDataTOM_importXML_axes.c
*/
int xDataXML_axesElememtToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_axes *axes );
int xDataXML_axesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_axes *axes );

/*
c Stuff in xDataTOM_importXML_XYs.c
*/
int xDataXML_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );
int xDataXML_XYsDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, int index, int length, double value, double accuracy,
        enum xDataTOM_subAxesType subAxesType, int axesOffest, xDataTOM_axes *axes, xDataTOM_interpolation *interpolation );
int xDataXML_stringToDoubles( statusMessageReporting *smr, xDataXML_element *XE, char const *s, int length, double *d );

/*
c Stuff in xDataTOM_importXML_regionsXYs.c
*/
int xDataXML_regionsXYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

/*
c Stuff in xDataTOM_importXML_W_XYs.c
*/
int xDataXML_W_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );
int xDataXML_W_XYsDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_W_XYs *W_XYs, int index, double value, xDataTOM_axes *axes,
    int axesOffset );

/*
c Stuff in xDataTOM_importXML_V_W_XYs.c
*/
int xDataXML_V_W_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

/*
c Stuff in xDataTOM_importXML_W_XYs_LegendreSeries.c
*/
int xDataXML_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );
int xDataXML_W_XYs_LegendreSeries_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_LegendreSeries *LegendreSeries );

/*
c Stuff in xDataTOM_importXML_regionsW_XYs_LegendreSeries.c
*/
int xDataXML_regionsW_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

/*
c Stuff in xDataTOM_importXML_V_W_XYs_LegendreSeries.c
*/
int xDataXML_V_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

/*
c Stuff in xDataTOM_importXML_polynomial.c
*/
int xDataXML_polynomialToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

/*
c Stuff in xDataTOM_importXML_KalbachMann.c
*/
int xDataXML_KalbachMannToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE );

#if defined __cplusplus
    }
    }
#endif

#endif              /* End of xDataTOM_importXML_private_h_included. */
