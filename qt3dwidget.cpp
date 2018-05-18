#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QCameraLens>
#include <Qt3DCore/QTransform>
#include <Qt3DCore/QAspectEngine>

#include <Qt3DInput/QInputAspect>

#include <Qt3DRender/QRenderAspect>
#include <Qt3DExtras/QForwardRenderer>

#include "qt3dwindow.h"
#include "qt3dwidget.h"
#include "viewer3d.h"



Qt3DWidget::Qt3DWidget(QWidget *parent,vector<Triangle>* meshVector)
   : QWidget(parent)
{
    Viewer3D*  m_viewer3D = new Viewer3D(this);
    m_viewer3D->sceneModifier()->addTriangleMeshCustomMaterial(QString(), meshVector);
    m_viewer3D->show();
}
