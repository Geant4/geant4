/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef MCGIDI_misc_h_included
#define MCGIDI_misc_h_included

#include <statusMessageReporting.h>
#include <xDataTOM_importXML_private.h>
#include "MCGIDI_private.h"

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

char const *MCGIDI_misc_pointerToTOMAttributeIfAllOk( statusMessageReporting *smr, char const *path, int required,
        xDataTOM_attributionList *attributes, char const *name, char const *file, int line );
char const *MCGIDI_misc_pointerToAttributeIfAllOk( statusMessageReporting *smr, xDataXML_element *element, char const *path, int required,
        xDataTOM_attributionList *attributes, char const *name, char const *file, int line );
int MCGIDI_misc_setMessageError_Element( statusMessageReporting *smr, void *userInterface, xDataXML_element *element, char const *file, int line, int code,
    char const *fmt, ... );
char *MCGIDI_misc_getAbsPath( statusMessageReporting *smr, char const *fileName );
int MCGIDI_misc_copyXMLAttributesToTOM( statusMessageReporting *smr, xDataTOM_attributionList *TOM, xDataXML_attributionList *XML );

#define MCGIDI_misc_pointerToTOMAttributeIfAllOk2( smr, required, attributes, name ) \
    MCGIDI_misc_pointerToTOMAttributeIfAllOk( smr, NULL, required, attributes, name, __FILE__, __LINE__ )
#define MCGIDI_misc_pointerToTOMAttributeIfAllOk3( smr, path, required, attributes, name ) \
    MCGIDI_misc_pointerToTOMAttributeIfAllOk( smr, path, required, attributes, name, __FILE__, __LINE__ )

#define MCGIDI_misc_pointerToAttributeIfAllOk2( smr, element, required, attributes, name ) \
    MCGIDI_misc_pointerToAttributeIfAllOk( smr, element, NULL, required, attributes, name, __FILE__, __LINE__ )
#define MCGIDI_misc_pointerToAttributeIfAllOk3( smr, path, required, attributes, name ) \
    MCGIDI_misc_pointerToAttributeIfAllOk( smr, NULL, path, required, attributes, name, __FILE__, __LINE__ )
enum xDataTOM_frame MCGIDI_misc_getProductFrame( statusMessageReporting *smr, xDataTOM_element *frameElement );

double MCGIDI_misc_getUnitConversionFactor( statusMessageReporting *smr, char const *fromUnit, char const *toUnit );
ptwXYPoints *MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( statusMessageReporting *smr, xDataTOM_XYs *XYs, 
        ptwXY_interpolation interpolation, char const *units[2] );
ptwXYPoints *MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( statusMessageReporting *smr, xDataTOM_element *linear, char const *toUnits[2] );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of MCGIDI_misc_h_included. */
