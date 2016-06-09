/****************************************************************************
** G4OpenGLQtViewer meta object code from reading C++ file 'G4OpenGLQtViewer.hh'
**
** Created: Mon Nov 12 10:34:24 2007
**      by: The Qt MOC ($Id: G4OpenGLQtViewer_moc.cc,v 1.5 2007/11/14 11:49:00 lgarnier Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#undef QT_NO_COMPAT
#include "../include/G4OpenGLQtViewer.hh"
#include <qmetaobject.h>
#include <qapplication.h>

#if QT_VERSION < 0x040000

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *G4OpenGLQtViewer::className() const
{
    return "G4OpenGLQtViewer";
}

QMetaObject *G4OpenGLQtViewer::metaObj = 0;
static QMetaObjectCleanUp cleanUp_G4OpenGLQtViewer( "G4OpenGLQtViewer", &G4OpenGLQtViewer::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString G4OpenGLQtViewer::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4OpenGLQtViewer", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString G4OpenGLQtViewer::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4OpenGLQtViewer", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* G4OpenGLQtViewer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QObject::staticMetaObject();
    static const QUMethod slot_0 = {"actionDrawingWireframe", 0, 0 };
    static const QUMethod slot_1 = {"actionDrawingLineRemoval", 0, 0 };
    static const QUMethod slot_2 = {"actionDrawingSurfaceRemoval", 0, 0 };
    static const QUMethod slot_3 = {"actionDrawingLineSurfaceRemoval", 0, 0 };
    static const QUMethod slot_4 = {"actionControlPanels", 0, 0 };
    static const QUMethod slot_5 = {"actionExitG4", 0, 0 };
    static const QUMethod slot_6 = {"actionCreateEPS", 0, 0 };
    static const QUParameter param_slot_7[] = {
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_7 = {"toggleDrawingAction", 1, param_slot_7 };
    static const QUParameter param_slot_8[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_8 = {"toggleMouseAction", 1, param_slot_8 };
    static const QUParameter param_slot_9[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_9 = {"toggleRepresentation", 1, param_slot_9 };
    static const QUParameter param_slot_10[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_10 = {"toggleBackground", 1, param_slot_10 };
    static const QUParameter param_slot_11[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_11 = {"toggleTransparency", 1, param_slot_11 };
    static const QUParameter param_slot_12[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_12 = {"toggleAntialiasing", 1, param_slot_12 };
    static const QUParameter param_slot_13[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_13 = {"toggleHaloing", 1, param_slot_13 };
    static const QUParameter param_slot_14[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_14 = {"toggleAux", 1, param_slot_14 };
    static const QUParameter param_slot_15[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_15 = {"toggleFullScreen", 1, param_slot_15 };
    static const QUMethod slot_16 = {"dialogClosed", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "actionDrawingWireframe()", &slot_0, QMetaData::Private },
	{ "actionDrawingLineRemoval()", &slot_1, QMetaData::Private },
	{ "actionDrawingSurfaceRemoval()", &slot_2, QMetaData::Private },
	{ "actionDrawingLineSurfaceRemoval()", &slot_3, QMetaData::Private },
	{ "actionControlPanels()", &slot_4, QMetaData::Private },
	{ "actionExitG4()", &slot_5, QMetaData::Private },
	{ "actionCreateEPS()", &slot_6, QMetaData::Private },
	{ "toggleDrawingAction(int)", &slot_7, QMetaData::Private },
	{ "toggleMouseAction(bool)", &slot_8, QMetaData::Private },
	{ "toggleRepresentation(bool)", &slot_9, QMetaData::Private },
	{ "toggleBackground(bool)", &slot_10, QMetaData::Private },
	{ "toggleTransparency(bool)", &slot_11, QMetaData::Private },
	{ "toggleAntialiasing(bool)", &slot_12, QMetaData::Private },
	{ "toggleHaloing(bool)", &slot_13, QMetaData::Private },
	{ "toggleAux(bool)", &slot_14, QMetaData::Private },
	{ "toggleFullScreen(bool)", &slot_15, QMetaData::Private },
	{ "dialogClosed()", &slot_16, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"G4OpenGLQtViewer", parentObject,
	slot_tbl, 17,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_G4OpenGLQtViewer.setMetaObject( metaObj );
    return metaObj;
}

void* G4OpenGLQtViewer::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "G4OpenGLQtViewer" ) )
	return this;
    if ( !qstrcmp( clname, "G4OpenGLViewer" ) )
	return (G4OpenGLViewer*)this;
    return QObject::qt_cast( clname );
}

bool G4OpenGLQtViewer::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: actionDrawingWireframe(); break;
    case 1: actionDrawingLineRemoval(); break;
    case 2: actionDrawingSurfaceRemoval(); break;
    case 3: actionDrawingLineSurfaceRemoval(); break;
    case 4: actionControlPanels(); break;
    case 5: actionExitG4(); break;
    case 6: actionCreateEPS(); break;
    case 7: toggleDrawingAction((int)static_QUType_int.get(_o+1)); break;
    case 8: toggleMouseAction((bool)static_QUType_bool.get(_o+1)); break;
    case 9: toggleRepresentation((bool)static_QUType_bool.get(_o+1)); break;
    case 10: toggleBackground((bool)static_QUType_bool.get(_o+1)); break;
    case 11: toggleTransparency((bool)static_QUType_bool.get(_o+1)); break;
    case 12: toggleAntialiasing((bool)static_QUType_bool.get(_o+1)); break;
    case 13: toggleHaloing((bool)static_QUType_bool.get(_o+1)); break;
    case 14: toggleAux((bool)static_QUType_bool.get(_o+1)); break;
    case 15: toggleFullScreen((bool)static_QUType_bool.get(_o+1)); break;
    case 16: dialogClosed(); break;
    default:
	return QObject::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool G4OpenGLQtViewer::qt_emit( int _id, QUObject* _o )
{
    return QObject::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool G4OpenGLQtViewer::qt_property( int id, int f, QVariant* v)
{
    return QObject::qt_property( id, f, v);
}

bool G4OpenGLQtViewer::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES

#endif

#endif


/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtViewer.hh'
**
** Created: Tue Sep 18 17:43:43 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if QT_VERSION >= 0x040000

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "../include/G4OpenGLQtViewer.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4OpenGLQtViewer.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.2.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_G4OpenGLQtViewer[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      17,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      17,   42,   42,   42, 0x08,
      43,   42,   42,   42, 0x08,
      70,   42,   42,   42, 0x08,
     100,   42,   42,   42, 0x08,
     134,   42,   42,   42, 0x08,
     156,   42,   42,   42, 0x08,
     171,   42,   42,   42, 0x08,
     189,   42,   42,   42, 0x08,
     214,   42,   42,   42, 0x08,
     238,   42,   42,   42, 0x08,
     265,   42,   42,   42, 0x08,
     288,   42,   42,   42, 0x08,
     313,   42,   42,   42, 0x08,
     338,   42,   42,   42, 0x08,
     358,   42,   42,   42, 0x08,
     374,   42,   42,   42, 0x08,
     397,   42,   42,   42, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_G4OpenGLQtViewer[] = {
    "G4OpenGLQtViewer\0actionDrawingWireframe()\0"
    "\0actionDrawingLineRemoval()\0"
    "actionDrawingSurfaceRemoval()\0"
    "actionDrawingLineSurfaceRemoval()\0"
    "actionControlPanels()\0actionExitG4()\0"
    "actionCreateEPS()\0toggleDrawingAction(int)\0"
    "toggleMouseAction(bool)\0"
    "toggleRepresentation(bool)\0"
    "toggleBackground(bool)\0toggleTransparency(bool)\0"
    "toggleAntialiasing(bool)\0toggleHaloing(bool)\0"
    "toggleAux(bool)\0toggleFullScreen(bool)\0"
    "dialogClosed()\0"
};

const QMetaObject G4OpenGLQtViewer::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_G4OpenGLQtViewer,
      qt_meta_data_G4OpenGLQtViewer, 0 }
};

const QMetaObject *G4OpenGLQtViewer::metaObject() const
{
    return &staticMetaObject;
}

void *G4OpenGLQtViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4OpenGLQtViewer))
	return static_cast<void*>(const_cast<G4OpenGLQtViewer*>(this));
    if (!strcmp(_clname, "G4OpenGLViewer"))
	return static_cast<G4OpenGLViewer*>(const_cast<G4OpenGLQtViewer*>(this));
    return QObject::qt_metacast(_clname);
}

int G4OpenGLQtViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: actionDrawingWireframe(); break;
        case 1: actionDrawingLineRemoval(); break;
        case 2: actionDrawingSurfaceRemoval(); break;
        case 3: actionDrawingLineSurfaceRemoval(); break;
        case 4: actionControlPanels(); break;
        case 5: actionExitG4(); break;
        case 6: actionCreateEPS(); break;
        case 7: toggleDrawingAction((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: toggleMouseAction((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: toggleRepresentation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: toggleBackground((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: toggleTransparency((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: toggleAntialiasing((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: toggleHaloing((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: toggleAux((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 15: toggleFullScreen((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: dialogClosed(); break;
        }
        _id -= 17;
    }
    return _id;
}
#endif
#endif
