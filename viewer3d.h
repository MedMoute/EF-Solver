#ifndef VIEWER3D_H
#define VIEWER3D_H

#include <QWidget>
#include <Qt3DCore/QEntity>
#include <QPointer>
#include <Qt3DExtras/Qt3DWindow>
#include <QDialog>

#include<QMatrix4x4>

#include "scenemodifier.h"

class Viewer3D : public QDialog
{
    Q_OBJECT

    public:
        Viewer3D(QWidget *parent = 0);
        SceneModifier* sceneModifier() {return m_sceneModifier;}

    private:

        QPointer<Qt3DCore::QEntity> m_rootEntity;
        QPointer<SceneModifier> m_sceneModifier;
        Qt3DExtras::Qt3DWindow *m_view;
        QPoint m_moveStartPoint;
        QMatrix4x4 m_cameraMatrix;
};
#endif // VIEWER3D_H
