/****************************************************************************
** G4UIQt meta object code from reading C++ file 'G4UIQt.hh'
**
** Created: Tue Nov 13 18:05:04 2007
**      by: The Qt MOC ($Id: G4UIQt_moc_030305.cc,v 1.1 2007-11-13 17:03:50 lgarnier Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#ifdef G4UI_BUILD_QT_SESSION

#undef QT_NO_COMPAT
#include "../include/G4UIQt.hh"
#include <qmetaobject.h>
#include <qapplication.h>

#if QT_VERSION < 0x040202

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *G4UIQt::className() const
{
    return "G4UIQt";
}

QMetaObject *G4UIQt::metaObj = 0;
static QMetaObjectCleanUp cleanUp_G4UIQt( "G4UIQt", &G4UIQt::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString G4UIQt::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4UIQt", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString G4UIQt::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "G4UIQt", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* G4UIQt::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QObject::staticMetaObject();
    static const QUMethod slot_0 = {"ClearButtonCallback", 0, 0 };
    static const QUMethod slot_1 = {"CommandEnteredCallback", 0, 0 };
    static const QUParameter param_slot_2[] = {
	{ 0, &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_2 = {"ButtonCallback", 1, param_slot_2 };
    static const QUMethod slot_3 = {"HelpTreeClicCallback", 0, 0 };
    static const QUParameter param_slot_4[] = {
	{ 0, &static_QUType_ptr, "QListViewItem", QUParameter::In },
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_4 = {"HelpTreeDoubleClicCallback", 2, param_slot_4 };
    static const QUParameter param_slot_5[] = {
	{ 0, &static_QUType_ptr, "QTreeWidgetItem", QUParameter::In },
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_5 = {"HelpTreeDoubleClicCallback", 2, param_slot_5 };
    static const QUMethod slot_6 = {"ShowHelpCallback", 0, 0 };
    static const QUMethod slot_7 = {"CommandHistoryCallback", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "ClearButtonCallback()", &slot_0, QMetaData::Private },
	{ "CommandEnteredCallback()", &slot_1, QMetaData::Private },
	{ "ButtonCallback(const QString&)", &slot_2, QMetaData::Private },
	{ "HelpTreeClicCallback()", &slot_3, QMetaData::Private },
	{ "HelpTreeDoubleClicCallback(QListViewItem*,int)", &slot_4, QMetaData::Private },
	{ "HelpTreeDoubleClicCallback(QTreeWidgetItem*,int)", &slot_5, QMetaData::Private },
	{ "ShowHelpCallback()", &slot_6, QMetaData::Private },
	{ "CommandHistoryCallback()", &slot_7, QMetaData::Private }
    };
    static const QUParameter param_signal_0[] = {
	{ "text", &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod signal_0 = {"myClicked", 1, param_signal_0 };
    static const QMetaData signal_tbl[] = {
	{ "myClicked(const QString&)", &signal_0, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"G4UIQt", parentObject,
	slot_tbl, 8,
	signal_tbl, 1,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_G4UIQt.setMetaObject( metaObj );
    return metaObj;
}

void* G4UIQt::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "G4UIQt" ) )
	return this;
    if ( !qstrcmp( clname, "G4VBasicShell" ) )
	return (G4VBasicShell*)this;
    if ( !qstrcmp( clname, "G4VInteractiveSession" ) )
	return (G4VInteractiveSession*)this;
    return QObject::qt_cast( clname );
}

// SIGNAL myClicked
void G4UIQt::myClicked( const QString& t0 )
{
    activate_signal( staticMetaObject()->signalOffset() + 0, t0 );
}

bool G4UIQt::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: ClearButtonCallback(); break;
    case 1: CommandEnteredCallback(); break;
    case 2: ButtonCallback((const QString&)static_QUType_QString.get(_o+1)); break;
    case 3: HelpTreeClicCallback(); break;
    case 4: HelpTreeDoubleClicCallback((QListViewItem*)static_QUType_ptr.get(_o+1),(int)static_QUType_int.get(_o+2)); break;
    case 5: HelpTreeDoubleClicCallback((QTreeWidgetItem*)static_QUType_ptr.get(_o+1),(int)static_QUType_int.get(_o+2)); break;
    case 6: ShowHelpCallback(); break;
    case 7: CommandHistoryCallback(); break;
    default:
	return QObject::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool G4UIQt::qt_emit( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->signalOffset() ) {
    case 0: myClicked((const QString&)static_QUType_QString.get(_o+1)); break;
    default:
	return QObject::qt_emit(_id,_o);
    }
    return TRUE;
}
#ifndef QT_NO_PROPERTIES

bool G4UIQt::qt_property( int id, int f, QVariant* v)
{
    return QObject::qt_property( id, f, v);
}

bool G4UIQt::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES

#endif
