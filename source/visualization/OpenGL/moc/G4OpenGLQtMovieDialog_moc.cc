/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtMovieDialog.hh'
**
** Created: Tue May 23 15:27:06 2017
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/G4OpenGLQtMovieDialog.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4OpenGLQtMovieDialog.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_G4OpenGLQtMovieDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      23,   22,   22,   22, 0x0a,
      41,   22,   22,   22, 0x0a,
      53,   22,   48,   22, 0x0a,
      80,   22,   48,   22, 0x0a,
     110,   22,   48,   22, 0x0a,
     138,   22,   22,   22, 0x08,
     164,   22,   22,   22, 0x08,
     187,   22,   22,   22, 0x08,
     214,   22,   22,   22, 0x08,
     231,   22,   22,   22, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_G4OpenGLQtMovieDialog[] = {
    "G4OpenGLQtMovieDialog\0\0stopFinishClose()\0"
    "save()\0bool\0checkEncoderSwParameters()\0"
    "checkSaveFileNameParameters()\0"
    "checkTempFolderParameters()\0"
    "selectEncoderPathAction()\0"
    "selectTempPathAction()\0"
    "selectSaveFileNameAction()\0resetRecording()\0"
    "enabledApplyButton()\0"
};

void G4OpenGLQtMovieDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        G4OpenGLQtMovieDialog *_t = static_cast<G4OpenGLQtMovieDialog *>(_o);
        switch (_id) {
        case 0: _t->stopFinishClose(); break;
        case 1: _t->save(); break;
        case 2: { bool _r = _t->checkEncoderSwParameters();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 3: { bool _r = _t->checkSaveFileNameParameters();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 4: { bool _r = _t->checkTempFolderParameters();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 5: _t->selectEncoderPathAction(); break;
        case 6: _t->selectTempPathAction(); break;
        case 7: _t->selectSaveFileNameAction(); break;
        case 8: _t->resetRecording(); break;
        case 9: _t->enabledApplyButton(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData G4OpenGLQtMovieDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject G4OpenGLQtMovieDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_G4OpenGLQtMovieDialog,
      qt_meta_data_G4OpenGLQtMovieDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &G4OpenGLQtMovieDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *G4OpenGLQtMovieDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *G4OpenGLQtMovieDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4OpenGLQtMovieDialog))
        return static_cast<void*>(const_cast< G4OpenGLQtMovieDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int G4OpenGLQtMovieDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
