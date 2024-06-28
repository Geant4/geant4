/********************************************************************************
** Form generated from reading UI file 'OIQtListsDialogwS1786.ui'
**
** Created by: Qt User Interface Compiler version 5.6.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef OIQTLISTSDIALOGWS1786_H
#define OIQTLISTSDIALOGWS1786_H

#include <QVariant>
#include <QAction>
#include <QApplication>
#include <QButtonGroup>
#include <QDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QWidget>

QT_BEGIN_NAMESPACE

class Ui_Dialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton_2;
    QPushButton *pushButton_3;
    QPushButton *pushButton;
    QWidget *layoutWidget;
    QFormLayout *formLayout;
    QLabel *label;
    QLineEdit *lineEdit;
    QListWidget *listWidget;
    QWidget *layoutWidget1;
    QFormLayout *formLayout_2;
    QLabel *label_2;
    QListWidget *listWidget1;

    void setupUi(QDialog *Dialog)
    {
        if (Dialog->objectName().isEmpty())
            Dialog->setObjectName(QStringLiteral("Dialog"));
        Dialog->resize(224, 554);
        QFont font;
        font.setPointSize(12);
        Dialog->setFont(font);
        horizontalLayoutWidget = new QWidget(Dialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(10, 260, 201, 26));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(8);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMaximumSize);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        pushButton_2 = new QPushButton(horizontalLayoutWidget);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pushButton_2->sizePolicy().hasHeightForWidth());
        pushButton_2->setSizePolicy(sizePolicy);
        pushButton_2->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));

        horizontalLayout->addWidget(pushButton_2);

        pushButton_3 = new QPushButton(horizontalLayoutWidget);
        pushButton_3->setObjectName(QStringLiteral("pushButton_3"));
        sizePolicy.setHeightForWidth(pushButton_3->sizePolicy().hasHeightForWidth());
        pushButton_3->setSizePolicy(sizePolicy);
        pushButton_3->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));

        horizontalLayout->addWidget(pushButton_3);

        pushButton = new QPushButton(horizontalLayoutWidget);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        sizePolicy.setHeightForWidth(pushButton->sizePolicy().hasHeightForWidth());
        pushButton->setSizePolicy(sizePolicy);
        pushButton->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));

        horizontalLayout->addWidget(pushButton);

        layoutWidget = new QWidget(Dialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 10, 201, 241));
        formLayout = new QFormLayout(layoutWidget);
        formLayout->setObjectName(QStringLiteral("formLayout"));
        formLayout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
        formLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(layoutWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));

        formLayout->setWidget(0, QFormLayout::LabelRole, label);

        lineEdit = new QLineEdit(layoutWidget);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));
        lineEdit->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));
        lineEdit->setReadOnly(true);

        formLayout->setWidget(1, QFormLayout::SpanningRole, lineEdit);

        listWidget = new QListWidget(layoutWidget);
        listWidget->setObjectName(QStringLiteral("listWidget"));
        listWidget->setFont(font);

        formLayout->setWidget(3, QFormLayout::SpanningRole, listWidget);

        layoutWidget1 = new QWidget(Dialog);
        layoutWidget1->setObjectName(QStringLiteral("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 300, 201, 241));
        formLayout_2 = new QFormLayout(layoutWidget1);
        formLayout_2->setObjectName(QStringLiteral("formLayout_2"));
        formLayout_2->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(layoutWidget1);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setStyleSheet(QStringLiteral("font: 12pt \"Sans Serif\";"));

        formLayout_2->setWidget(0, QFormLayout::LabelRole, label_2);

        listWidget1 = new QListWidget(layoutWidget1);
        listWidget1->setObjectName(QStringLiteral("listWidget1"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(1);
        sizePolicy1.setVerticalStretch(1);
        sizePolicy1.setHeightForWidth(listWidget1->sizePolicy().hasHeightForWidth());
        listWidget1->setSizePolicy(sizePolicy1);
        listWidget1->setFont(font);

        formLayout_2->setWidget(1, QFormLayout::SpanningRole, listWidget1);


        retranslateUi(Dialog);

        QMetaObject::connectSlotsByName(Dialog);
    } // setupUi

    void retranslateUi(QDialog *Dialog)
    {
        Dialog->setWindowTitle(QApplication::translate("Dialog", "OIQt", 0));
#ifndef QT_NO_TOOLTIP
        pushButton_2->setToolTip(QApplication::translate("Dialog", "Delete bookmark", 0));
#endif // QT_NO_TOOLTIP
        pushButton_2->setText(QApplication::translate("Dialog", "Del", 0));
#ifndef QT_NO_TOOLTIP
        pushButton_3->setToolTip(QApplication::translate("Dialog", "Rename bookmark", 0));
#endif // QT_NO_TOOLTIP
        pushButton_3->setText(QApplication::translate("Dialog", "Ren", 0));
#ifndef QT_NO_TOOLTIP
        pushButton->setToolTip(QApplication::translate("Dialog", "Sort bookmarks", 0));
#endif // QT_NO_TOOLTIP
        pushButton->setText(QApplication::translate("Dialog", "Sort", 0));
        label->setText(QApplication::translate("Dialog", "Bookmarks", 0));
        label_2->setText(QApplication::translate("Dialog", "Elements   [S mm]", 0));
    } // retranslateUi

};

namespace Ui {
    class Dialog: public Ui_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // OIQTLISTSDIALOGWS1786_H
