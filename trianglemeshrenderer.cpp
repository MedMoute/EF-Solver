#include "trianglemeshrenderer.h"

TriangleMeshRenderer::TriangleMeshRenderer(const std::vector<Triangle> *meshVector, QNode *parent)
    : Qt3DRender::QGeometryRenderer(parent)
{
    setPrimitiveType(Qt3DRender::QGeometryRenderer::Triangles);
    setGeometry(new TriangleMeshGeometry(meshVector, this));
}

TriangleMeshRenderer::~TriangleMeshRenderer()
{
}

TriangleMeshGeometry::TriangleMeshGeometry(const std::vector<Triangle>* meshVector, TriangleMeshRenderer *parent)
    : Qt3DRender::QGeometry(parent)
{
    Qt3DRender::QBuffer *vertexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::VertexBuffer, this);
    Qt3DRender::QBuffer *indexDataBuffer = new Qt3DRender::QBuffer(Qt3DRender::QBuffer::IndexBuffer, this);

    // Vertexbuffer
    QByteArray vertexBufferData;
    // Buffer size = triangle count * 3 * (3 + 3 + 3), 3 vertices per trinalge, each 3 floats for vertex position x,y,z, 3 floats normal and 3 floats color
    int bytesPerVertex = 9 * sizeof(float);
    int bytesPerTriangle = 3 * bytesPerVertex;
    vertexBufferData.resize(static_cast<int>(meshVector->size()) * bytesPerTriangle);
    char* pByte = vertexBufferData.data();
    int i = 0;
    // Indexbuffer
    QByteArray indexBufferData;
    indexBufferData.resize(static_cast<int>(meshVector->size()) * 3 * sizeof(uint));
    uint* rawIndexArray = reinterpret_cast<uint*>(indexBufferData.data());
    int idx = 0;

    for (int n = 0; n < meshVector->size(); ++n)
    {
        QVector3D nt = QVector3D::normal((*meshVector)[n].vertices[0].p, (*meshVector)[n].vertices[1].p, (*meshVector)[n].vertices[2].p);

        for (int v = 0; v < 3; ++v)
        {
            // Vertex
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].p.x(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].p.y(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].p.z(); pByte += 4;
            // Normal
            *reinterpret_cast<float*>(pByte) = nt.x(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = nt.y(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = nt.z(); pByte += 4;
            // Color
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].c.x(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].c.y(); pByte += 4;
            *reinterpret_cast<float*>(pByte) = (*meshVector)[n].vertices[v].c.z(); pByte += 4;

            // Index
            rawIndexArray[idx] = static_cast<uint>(idx++);
        }
    }

    vertexDataBuffer->setData(vertexBufferData);
    indexDataBuffer->setData(indexBufferData);

    // Attributes
    Qt3DRender::QAttribute *positionAttribute = new Qt3DRender::QAttribute();
    positionAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    positionAttribute->setBuffer(vertexDataBuffer);
    positionAttribute->setDataType(Qt3DRender::QAttribute::Float);
    positionAttribute->setDataSize(3);
    positionAttribute->setByteOffset(0);
    positionAttribute->setByteStride(bytesPerVertex);
    positionAttribute->setCount(3 * static_cast<int>(meshVector->size()));
    positionAttribute->setName(Qt3DRender::QAttribute::defaultPositionAttributeName());

    Qt3DRender::QAttribute *normalAttribute = new Qt3DRender::QAttribute();
    normalAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    normalAttribute->setBuffer(vertexDataBuffer);
    normalAttribute->setDataType(Qt3DRender::QAttribute::Float);
    normalAttribute->setDataSize(3);
    normalAttribute->setByteOffset(3 * sizeof(float));
    normalAttribute->setByteStride(bytesPerVertex);
    normalAttribute->setCount(3 * static_cast<int>(meshVector->size()));
    normalAttribute->setName(Qt3DRender::QAttribute::defaultNormalAttributeName());

    Qt3DRender::QAttribute *colorAttribute = new Qt3DRender::QAttribute();
    colorAttribute->setAttributeType(Qt3DRender::QAttribute::VertexAttribute);
    colorAttribute->setBuffer(vertexDataBuffer);
    colorAttribute->setDataType(Qt3DRender::QAttribute::Float);
    colorAttribute->setDataSize(3);
    colorAttribute->setByteOffset(6 * sizeof(float));
    colorAttribute->setByteStride(bytesPerVertex);
    colorAttribute->setCount(3 * static_cast<int>(meshVector->size()));
    colorAttribute->setName(Qt3DRender::QAttribute::defaultColorAttributeName());

    Qt3DRender::QAttribute *indexAttribute = new Qt3DRender::QAttribute();
    indexAttribute->setAttributeType(Qt3DRender::QAttribute::IndexAttribute);
    indexAttribute->setBuffer(indexDataBuffer);
    indexAttribute->setDataType(Qt3DRender::QAttribute::UnsignedInt);
    indexAttribute->setDataSize(1);
    indexAttribute->setByteOffset(0);
    indexAttribute->setByteStride(0);
    indexAttribute->setCount(3 * static_cast<int>(meshVector->size()));

    addAttribute(positionAttribute);
    addAttribute(normalAttribute);
    addAttribute(colorAttribute);
    addAttribute(indexAttribute);

    parent->setGeometry(this);
}
