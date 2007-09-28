/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtViewer.hh'
**
** Created: Tue Sep 18 17:43:43 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

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
