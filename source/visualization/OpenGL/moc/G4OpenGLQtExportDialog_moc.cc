/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtExportDialog.hh'
**
** Created: Tue May 23 15:27:06 2017
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/G4OpenGLQtExportDialog.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'G4OpenGLQtExportDialog.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_G4OpenGLQtExportDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      24,   23,   23,   23, 0x0a,
      40,   23,   23,   23, 0x0a,
      58,   23,   23,   23, 0x0a,
      84,   23,   23,   23, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_G4OpenGLQtExportDialog[] = {
    "G4OpenGLQtExportDialog\0\0changeSizeBox()\0"
    "changeVectorEPS()\0textWidthChanged(QString)\0"
    "textHeightChanged(QString)\0"
};

void G4OpenGLQtExportDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        G4OpenGLQtExportDialog *_t = static_cast<G4OpenGLQtExportDialog *>(_o);
        switch (_id) {
        case 0: _t->changeSizeBox(); break;
        case 1: _t->changeVectorEPS(); break;
        case 2: _t->textWidthChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: _t->textHeightChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData G4OpenGLQtExportDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject G4OpenGLQtExportDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_G4OpenGLQtExportDialog,
      qt_meta_data_G4OpenGLQtExportDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &G4OpenGLQtExportDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *G4OpenGLQtExportDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *G4OpenGLQtExportDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_G4OpenGLQtExportDialog))
        return static_cast<void*>(const_cast< G4OpenGLQtExportDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int G4OpenGLQtExportDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
