#ifndef SCENEMODIFIER_H
#define SCENEMODIFIER_H
#include <QObject>
#include <Qt3DCore/QEntity>
#include "utils.h"

class SceneModifier : public QObject
{
    Q_OBJECT

    public:
        SceneModifier(Qt3DCore::QEntity* rootEntity);
        void addTriangleMeshCustomMaterial(QString name, const std::vector<Triangle> *meshVector);

    private:
        Qt3DCore::QEntity* m_rootEntity;
};
#endif // SCENEMODIFIER_H
