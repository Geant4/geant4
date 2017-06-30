/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtViewer.hh'
**
** Created: Tue May 23 15:27:06 2017
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/G4OpenGLQtViewer.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4OpenGLQtViewer.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_G4OpenGLQtViewer[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      29,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   17,   17,   17, 0x0a,
      36,   17,   17,   17, 0x09,
      71,   17,   17,   17, 0x08,
      89,   17,   17,   17, 0x08,
     119,   17,   17,   17, 0x08,
     143,   17,   17,   17, 0x08,
     170,   17,   17,   17, 0x08,
     194,   17,   17,   17, 0x08,
     210,   17,   17,   17, 0x08,
     233,   17,   17,   17, 0x08,
     258,   17,   17,   17, 0x08,
     281,   17,   17,   17, 0x08,
     306,   17,   17,   17, 0x08,
     331,   17,   17,   17, 0x08,
     351,   17,   17,   17, 0x08,
     367,   17,   17,   17, 0x08,
     393,   17,   17,   17, 0x08,
     416,   17,   17,   17, 0x08,
     440,   17,   17,   17, 0x08,
     465,   17,   17,   17, 0x08,
     495,  487,   17,   17, 0x08,
     547,   17,   17,   17, 0x08,
     588,   17,   17,   17, 0x08,
     604,   17,   17,   17, 0x08,
     629,   17,   17,   17, 0x08,
     658,   17,   17,   17, 0x08,
     686,   17,   17,   17, 0x08,
     719,  710,   17,   17, 0x08,
     768,   17,   17,   17, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_G4OpenGLQtViewer[] = {
    "G4OpenGLQtViewer\0\0startPauseVideo()\0"
    "updateToolbarAndMouseContextMenu()\0"
    "actionSaveImage()\0actionChangeBackgroundColor()\0"
    "actionChangeTextColor()\0"
    "actionChangeDefaultColor()\0"
    "actionMovieParameters()\0showShortcuts()\0"
    "toggleMouseAction(int)\0toggleSurfaceAction(int)\0"
    "toggleProjection(bool)\0toggleTransparency(bool)\0"
    "toggleAntialiasing(bool)\0toggleHaloing(bool)\0"
    "toggleAux(bool)\0toggleHiddenMarkers(bool)\0"
    "toggleFullScreen(bool)\0processEncodeFinished()\0"
    "processLookForFinished()\0processEncodeStdout()\0"
    "item,id\0sceneTreeComponentItemChanged(QTreeWidgetItem*,int)\0"
    "toggleSceneTreeComponentPickingCout(int)\0"
    "togglePicking()\0currentTabActivated(int)\0"
    "sceneTreeComponentSelected()\0"
    "changeDepthInSceneTree(int)\0"
    "changeSearchSelection()\0item,val\0"
    "changeColorAndTransparency(QTreeWidgetItem*,int)\0"
    "tableWidgetViewerSetItemChanged(QTableWidgetItem*)\0"
};

void G4OpenGLQtViewer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        G4OpenGLQtViewer *_t = static_cast<G4OpenGLQtViewer *>(_o);
        switch (_id) {
        case 0: _t->startPauseVideo(); break;
        case 1: _t->updateToolbarAndMouseContextMenu(); break;
        case 2: _t->actionSaveImage(); break;
        case 3: _t->actionChangeBackgroundColor(); break;
        case 4: _t->actionChangeTextColor(); break;
        case 5: _t->actionChangeDefaultColor(); break;
        case 6: _t->actionMovieParameters(); break;
        case 7: _t->showShortcuts(); break;
        case 8: _t->toggleMouseAction((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->toggleSurfaceAction((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->toggleProjection((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: _t->toggleTransparency((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: _t->toggleAntialiasing((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->toggleHaloing((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->toggleAux((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 15: _t->toggleHiddenMarkers((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: _t->toggleFullScreen((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->processEncodeFinished(); break;
        case 18: _t->processLookForFinished(); break;
        case 19: _t->processEncodeStdout(); break;
        case 20: _t->sceneTreeComponentItemChanged((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 21: _t->toggleSceneTreeComponentPickingCout((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 22: _t->togglePicking(); break;
        case 23: _t->currentTabActivated((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 24: _t->sceneTreeComponentSelected(); break;
        case 25: _t->changeDepthInSceneTree((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 26: _t->changeSearchSelection(); break;
        case 27: _t->changeColorAndTransparency((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 28: _t->tableWidgetViewerSetItemChanged((*reinterpret_cast< QTableWidgetItem*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData G4OpenGLQtViewer::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject G4OpenGLQtViewer::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_G4OpenGLQtViewer,
      qt_meta_data_G4OpenGLQtViewer, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &G4OpenGLQtViewer::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *G4OpenGLQtViewer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *G4OpenGLQtViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4OpenGLQtViewer))
        return static_cast<void*>(const_cast< G4OpenGLQtViewer*>(this));
    if (!strcmp(_clname, "G4OpenGLViewer"))
        return static_cast< G4OpenGLViewer*>(const_cast< G4OpenGLQtViewer*>(this));
    return QObject::qt_metacast(_clname);
}

int G4OpenGLQtViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 29)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 29;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
