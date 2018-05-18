#ifndef QT3DWIDGET_H
#define QT3DWIDGET_H

#include <QObject>
#include <QWidget>
#include <Qt3DExtras/Qt3DWindow>
#include "utils.h"

class Qt3DWidget
      : public QWidget
{
    Q_OBJECT
    QWidget *container;

public:
    explicit Qt3DWidget(QWidget *parent = nullptr,std::vector<Triangle>* tris);

};

#endif // QT3DWIDGET_H
