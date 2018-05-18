#ifndef TRIANGLEMESHRENDERER_H
#define TRIANGLEMESHRENDERER_H
#include<Qt3DRender>
#include "utils.h"

class TriangleMeshRenderer : public Qt3DRender::QGeometryRenderer
{
    Q_OBJECT
public:
    explicit TriangleMeshRenderer(const std::vector<Triangle>* meshVector, Qt3DCore::QNode *parent = 0);
    ~TriangleMeshRenderer();
};


class TriangleMeshGeometry : public Qt3DRender::QGeometry
{
    Q_OBJECT
public:
    TriangleMeshGeometry(const std::vector<Triangle>* meshVector, TriangleMeshRenderer *parent);
};


#endif // TRIANGLEMESHRENDERER_H
