/****************************************************************************
** Meta object code from reading C++ file 'G4UIQt.hh'
**
** Created: Thu Apr 20 12:02:10 2017
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/G4UIQt.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4UIQt.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_G4UIQt[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      25,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
       8,    7,    7,    7, 0x08,
      22,    7,    7,    7, 0x08,
      44,    7,    7,    7, 0x08,
      65,    7,    7,    7, 0x08,
      95,   90,    7,    7, 0x08,
     126,    7,    7,    7, 0x08,
     150,    7,    7,    7, 0x08,
     173,    7,    7,    7, 0x08,
     202,    7,    7,    7, 0x08,
     221,    7,    7,    7, 0x08,
     246,    7,    7,    7, 0x08,
     274,    7,    7,    7, 0x08,
     295,    7,    7,    7, 0x08,
     326,    7,    7,    7, 0x08,
     354,    7,    7,    7, 0x08,
     382,    7,    7,    7, 0x08,
     404,    7,    7,    7, 0x08,
     426,    7,    7,    7, 0x08,
     457,    7,    7,    7, 0x08,
     487,    7,    7,    7, 0x08,
     515,    7,    7,    7, 0x08,
     543,    7,    7,    7, 0x08,
     569,    7,    7,    7, 0x08,
     595,    7,    7,    7, 0x08,
     629,    7,    7,    7, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_G4UIQt[] = {
    "G4UIQt\0\0ExitSession()\0ClearButtonCallback()\0"
    "SaveOutputCallback()\0CommandEnteredCallback()\0"
    "text\0CommandEditedCallback(QString)\0"
    "ButtonCallback(QString)\0HelpTreeClicCallback()\0"
    "HelpTreeDoubleClicCallback()\0"
    "ShowHelpCallback()\0CommandHistoryCallback()\0"
    "LookForHelpStringCallback()\0"
    "UpdateTabWidget(int)\0"
    "ResizeTabWidget(QResizeEvent*)\0"
    "CoutFilterCallback(QString)\0"
    "ThreadComboBoxCallback(int)\0"
    "TabCloseCallback(int)\0ToolBoxActivated(int)\0"
    "VisParameterCallback(QWidget*)\0"
    "ChangeColorCallback(QWidget*)\0"
    "ChangeCursorAction(QString)\0"
    "ChangeSurfaceStyle(QString)\0"
    "OpenIconCallback(QString)\0"
    "SaveIconCallback(QString)\0"
    "ViewerPropertiesIconCallback(int)\0"
    "ChangePerspectiveOrtho(QString)\0"
};

void G4UIQt::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        G4UIQt *_t = static_cast<G4UIQt *>(_o);
        switch (_id) {
        case 0: _t->ExitSession(); break;
        case 1: _t->ClearButtonCallback(); break;
        case 2: _t->SaveOutputCallback(); break;
        case 3: _t->CommandEnteredCallback(); break;
        case 4: _t->CommandEditedCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->ButtonCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: _t->HelpTreeClicCallback(); break;
        case 7: _t->HelpTreeDoubleClicCallback(); break;
        case 8: _t->ShowHelpCallback(); break;
        case 9: _t->CommandHistoryCallback(); break;
        case 10: _t->LookForHelpStringCallback(); break;
        case 11: _t->UpdateTabWidget((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->ResizeTabWidget((*reinterpret_cast< QResizeEvent*(*)>(_a[1]))); break;
        case 13: _t->CoutFilterCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: _t->ThreadComboBoxCallback((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: _t->TabCloseCallback((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 16: _t->ToolBoxActivated((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: _t->VisParameterCallback((*reinterpret_cast< QWidget*(*)>(_a[1]))); break;
        case 18: _t->ChangeColorCallback((*reinterpret_cast< QWidget*(*)>(_a[1]))); break;
        case 19: _t->ChangeCursorAction((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 20: _t->ChangeSurfaceStyle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 21: _t->OpenIconCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 22: _t->SaveIconCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 23: _t->ViewerPropertiesIconCallback((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 24: _t->ChangePerspectiveOrtho((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData G4UIQt::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject G4UIQt::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_G4UIQt,
      qt_meta_data_G4UIQt, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &G4UIQt::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *G4UIQt::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *G4UIQt::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4UIQt))
        return static_cast<void*>(const_cast< G4UIQt*>(this));
    if (!strcmp(_clname, "G4VBasicShell"))
        return static_cast< G4VBasicShell*>(const_cast< G4UIQt*>(this));
    if (!strcmp(_clname, "G4VInteractiveSession"))
        return static_cast< G4VInteractiveSession*>(const_cast< G4UIQt*>(this));
    return QObject::qt_metacast(_clname);
}

int G4UIQt::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 25)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 25;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
