#if QT_VERSION > 0x040000

/****************************************************************************
** Meta object code from reading C++ file 'G4OpenGLQtExportDialog.hh'
**
** Created: Fri Sep 28 12:35:41 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

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
