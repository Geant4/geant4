/****************************************************************************
** Meta object code from reading C++ file 'G4UIQt.hh'
**
** Created: Mon Sep 17 14:28:15 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#ifdef G4UI_BUILD_QT_SESSION

#include "../include/G4UIQt.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4UIQt.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.2.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_G4UIQt[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
       7,   26,   31,   31, 0x05,

 // slots: signature, parameters, type, tag, flags
      32,   31,   31,   31, 0x08,
      54,   31,   31,   31, 0x08,
      79,   31,   31,   31, 0x08,
     103,   31,   31,   31, 0x08,
     126,  175,   31,   31, 0x08,
     177,   31,   31,   31, 0x08,
     196,   31,   31,   31, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_G4UIQt[] = {
    "G4UIQt\0myClicked(QString)\0text\0\0"
    "ClearButtonCallback()\0CommandEnteredCallback()\0"
    "ButtonCallback(QString)\0HelpTreeClicCallback()\0"
    "HelpTreeDoubleClicCallback(QTreeWidgetItem*,int)\0"
    ",\0ShowHelpCallback()\0CommandHistoryCallback()\0"
};

const QMetaObject G4UIQt::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_G4UIQt,
      qt_meta_data_G4UIQt, 0 }
};

const QMetaObject *G4UIQt::metaObject() const
{
    return &staticMetaObject;
}

void *G4UIQt::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4UIQt))
	return static_cast<void*>(const_cast<G4UIQt*>(this));
    if (!strcmp(_clname, "G4VBasicShell"))
	return static_cast<G4VBasicShell*>(const_cast<G4UIQt*>(this));
    if (!strcmp(_clname, "G4VInteractiveSession"))
	return static_cast<G4VInteractiveSession*>(const_cast<G4UIQt*>(this));
    return QObject::qt_metacast(_clname);
}

int G4UIQt::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: myClicked((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: ClearButtonCallback(); break;
        case 2: CommandEnteredCallback(); break;
        case 3: ButtonCallback((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: HelpTreeClicCallback(); break;
        case 5: HelpTreeDoubleClicCallback((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 6: ShowHelpCallback(); break;
        case 7: CommandHistoryCallback(); break;
        }
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void G4UIQt::myClicked(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

#endif
