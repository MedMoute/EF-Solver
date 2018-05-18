#include "viewer3d.h"
#include <Qt3DExtras/QForwardRenderer>
#include <QScreen>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <Qt3DRender/QCamera>
#include <Qt3DCore/QTransform>
#include <QWheelEvent>
#include <QEvent>

#include <Qt3DExtras/QOrbitCameraController>
#include "orbittransformcontroller.h"

Viewer3D::Viewer3D(QWidget *parent) :
    QDialog(parent)
{
    setAttribute(Qt::WA_DeleteOnClose);
    m_moveStartPoint.setX(-1);

    m_view = new Qt3DExtras::Qt3DWindow();

    m_view->installEventFilter(this);

    m_view->defaultFrameGraph()->setClearColor(QColor(QRgb(0x4d4d4f)));

    QWidget *container = QWidget::createWindowContainer(m_view);
    QSize screenSize = m_view->screen()->size();
    container->setMinimumSize(QSize(200, 100));
    container->setMaximumSize(screenSize);

    QHBoxLayout *hLayout = new QHBoxLayout(this);
    QVBoxLayout *vLayout = new QVBoxLayout();
    hLayout->addWidget(container, 1);

    setWindowTitle(QStringLiteral("Mesh Viewer"));

    // Root entity
    m_rootEntity = new Qt3DCore::QEntity();

    // Scene modifier
    m_sceneModifier = new SceneModifier(m_rootEntity);

    // Window geometry
    resize(parent->geometry().width() * 0.8, parent->geometry().height() * 0.8);
    move(parent->geometry().center() - QPoint(width() / 2, height() / 2));

    // Camera
    Qt3DRender::QCamera *camera = m_view->camera();
    camera->lens()->setPerspectiveProjection(45.0f, 16.0f/9.0f, 0.1f, 1000.0f);
    camera->setPosition(QVector3D(0, 0, 40.0f));
    camera->setViewCenter(QVector3D(0, 0, 0));

    // For camera controls
    Qt3DExtras::QOrbitCameraController *camController = new Qt3DExtras::QOrbitCameraController(m_rootEntity);
    camController->setLinearSpeed( 50.0f );
    camController->setLookSpeed( 180.0f );
    camController->setCamera(camera);

    // Set root object of the scene
    m_view->setRootEntity(m_rootEntity);
}

