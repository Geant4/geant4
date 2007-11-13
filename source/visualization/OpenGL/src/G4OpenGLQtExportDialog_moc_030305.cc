/****************************************************************************
** G4OpenGLQtExportDialog meta object code from reading C++ file 'G4OpenGLQtExportDialog.hh'
**
** Created: Mon Nov 12 10:36:59 2007
**      by: The Qt MOC ($Id: G4OpenGLQtExportDialog_moc_030305.cc,v 1.2 2007-11-13 15:56:11 lgarnier Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#undef QT_NO_COMPAT
#include "../include/G4OpenGLQtExportDialog.hh"
#include <qmetaobject.h>
#include <qapplication.h>

#if QT_VERSION < 0x040202

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
