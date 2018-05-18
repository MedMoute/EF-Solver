#ifndef UTILS_H
#define UTILS_H

#include<cmath>
#include<algorithm>
#include<tuple>
#include <QVector3D>

#include "R2.h"
#include "R3.h"

using namespace std;

struct Point
{
    QVector3D p; //point x, y, z
    QVector3D c; //color red, green, blue

    Point() {}

    Point(float xp, float yp, float zp)
    {
        p = QVector3D(xp, yp, zp);
        c = QVector3D(0, 0, 0);
    }
    Point(QVector3D pos, unsigned char r, unsigned char g, unsigned char b)
    {
        p = pos;
        c = QVector3D(static_cast<float>(r) / 255.f,
                      static_cast<float>(g) / 255.f,
                      static_cast<float>(b) / 255.f);
    }
};

struct Triangle
{
    Point vertices[3];

    Triangle()
    {
    }

    Triangle(Point p1, Point p2, Point p3)
    {
        vertices[0] = p1;
        vertices[1] = p2;
        vertices[2] = p3;
    }

};

class util
{
public:
   static float euclDistance( R2& a , R2& b);
   static float euclDistance( R3& a , R3& b);

   static float scalPdt( R2& a,  R2& b);
   static float scalPdt( R3& a,  R3& b);

   static int type2nnodes(int type);

   static void print_separator();

   static int factorial(int n);
   static int kron(int i,int j);

   static double simplex3Mesure(R3 a,R3 b,R3 c,R3 d);
   static double simplex2Mesure(R3 a,R3 b,R3 c);

};

#endif // UTILS_H
