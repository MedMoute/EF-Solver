#include "SceneModifier.h"
#include "TriangleMeshRenderer.h"
#include <Qt3DExtras/QPhongMaterial>

SceneModifier::SceneModifier(Qt3DCore::QEntity* rootEntity) :
    m_rootEntity(rootEntity),
    QObject(rootEntity)
{
}

void SceneModifier::addTriangleMeshCustomMaterial(QString name, const std::vector<Triangle>* meshVector)
{
    if (!m_rootEntity)
    {
        return;
    }

    // Mesh entity
    Qt3DCore::QEntity *triangleMeshEntity = new Qt3DCore::QEntity(m_rootEntity);
    triangleMeshEntity->setObjectName(QStringLiteral("customMeshEntity"));

    TriangleMeshRenderer *triangleMeshRenderer = new TriangleMeshRenderer(meshVector);
    Qt3DRender::QMaterial *material = new Qt3DExtras::QPhongMaterial(triangleMeshEntity);
    Qt3DCore::QTransform *transform = new Qt3DCore::QTransform;
    transform->setScale(1.f);

    triangleMeshEntity->addComponent(triangleMeshRenderer);
    triangleMeshEntity->addComponent(transform);
    triangleMeshEntity->addComponent(material);

    //emit meshAdded(name, triangleMeshEntity);
}
