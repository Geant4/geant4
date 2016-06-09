/****************************************************************************
** G4OpenGLQtExportDialog meta object code from reading C++ file 'G4OpenGLQtExportDialog.hh'
**
** Created: Mon Nov 12 10:36:59 2007
**      by: The Qt MOC ($Id: G4OpenGLQtExportDialog_moc.cc,v 1.5 2007/11/14 11:49:00 lgarnier Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#undef QT_NO_COMPAT
#include "../include/G4OpenGLQtExportDialog.hh"
#include <qmetaobject.h>
#include <qapplication.h>

#if QT_VERSION < 0x040000

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *G4OpenGLQtExportDialog::className() const
{
    return "G4OpenGLQtExportDialog";
}

QMetaObject *G4OpenGLQtExportDialog::metaObj = 0;
static QMetaObjectCleanUp cleanUp_G4OpenGLQtExportDialog( "G4OpenGLQtExportDialog", &G4OpenGLQtExportDialog::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString G4OpenGLQtExportDialog::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4OpenGLQtExportDialog", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString G4OpenGLQtExportDialog::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4OpenGLQtExportDialog", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* G4OpenGLQtExportDialog::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QDialog::staticMetaObject();
    static const QUParameter param_slot_0[] = {
	{ 0, &static_QUType_bool, 0, QUParameter::In }
    };
    static const QUMethod slot_0 = {"changeSizeBox", 1, param_slot_0 };
    static const QUParameter param_slot_1[] = {
	{ 0, &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_1 = {"textWidthChanged", 1, param_slot_1 };
    static const QUParameter param_slot_2[] = {
	{ 0, &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_2 = {"textHeightChanged", 1, param_slot_2 };
    static const QMetaData slot_tbl[] = {
	{ "changeSizeBox(bool)", &slot_0, QMetaData::Public },
	{ "textWidthChanged(const QString&)", &slot_1, QMetaData::Public },
	{ "textHeightChanged(const QString&)", &slot_2, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"G4OpenGLQtExportDialog", parentObject,
	slot_tbl, 3,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_G4OpenGLQtExportDialog.setMetaObject( metaObj );
    return metaObj;
}

void* G4OpenGLQtExportDialog::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "G4OpenGLQtExportDialog" ) )
	return this;
    return QDialog::qt_cast( clname );
}

bool G4OpenGLQtExportDialog::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: changeSizeBox((bool)static_QUType_bool.get(_o+1)); break;
    case 1: textWidthChanged((const QString&)static_QUType_QString.get(_o+1)); break;
    case 2: textHeightChanged((const QString&)static_QUType_QString.get(_o+1)); break;
    default:
	return QDialog::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool G4OpenGLQtExportDialog::qt_emit( int _id, QUObject* _o )
{
    return QDialog::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool G4OpenGLQtExportDialog::qt_property( int id, int f, QVariant* v)
{
    return QDialog::qt_property( id, f, v);
}

bool G4OpenGLQtExportDialog::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES

#endif
#endif


/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtExportDialog.hh'
**
** Created: Fri Sep 28 12:35:41 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if QT_VERSION >= 0x040000

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "../include/G4OpenGLQtExportDialog.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4OpenGLQtExportDialog.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.2.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_G4OpenGLQtExportDialog[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      23,   43,   43,   43, 0x0a,
      44,   43,   43,   43, 0x0a,
      70,   43,   43,   43, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_G4OpenGLQtExportDialog[] = {
    "G4OpenGLQtExportDialog\0changeSizeBox(bool)\0"
    "\0textWidthChanged(QString)\0"
    "textHeightChanged(QString)\0"
};

const QMetaObject G4OpenGLQtExportDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_G4OpenGLQtExportDialog,
      qt_meta_data_G4OpenGLQtExportDialog, 0 }
};

const QMetaObject *G4OpenGLQtExportDialog::metaObject() const
{
    return &staticMetaObject;
}

void *G4OpenGLQtExportDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4OpenGLQtExportDialog))
	return static_cast<void*>(const_cast<G4OpenGLQtExportDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int G4OpenGLQtExportDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: changeSizeBox((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: textWidthChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: textHeightChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 3;
    }
    return _id;
}
#endif
#endif
