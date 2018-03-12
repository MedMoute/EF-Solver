#include "r3.h"

R3::R3()
{
    x=0;
    y=0;
    z=0;
}

R3::R3(float _x, float _y,float _z)
{
    x=_x;
    y=_y;
    z=_z;

}

R3::R3(const R3& point)
{
    std::tuple<float,float,float> data=point.getCoords();
    x=std::get<0>(data);
    y=std::get<1>(data);
    z=std::get<2>(data);
}

std::tuple<float,float,float> R3::getCoords() const
{
    return std::make_tuple(x,y,z);
}

float R3::getX() const
{
    return x;
}

float R3::getY() const
{
    return y;
}

float R3::getZ() const
{
    return z;
}

void R3::setX(float _x)
{
    x=_x;
}

void R3::setY(float _y)
{
    y=_y;
}

void R3::setZ(float _z)
{
    z=_z;
}

void R3::displayCoords()
{
    std::cout<< "Point data : x ="<<getX()<<" --- y = "<<getY()<<" --- z = "<<getZ()<<std::endl;
}

bool R3::operator ==(const R3& b)
{
    return ((b.getX()==x)&&(b.getY()==y)&&b.getZ()==z);
}

std::ostream& R3::operator<<(std::ostream& out)
{
    out << "(" << x << ", " << y << ", " << z << ")";
    return out;
}

void R3::operator =(const R3& b)
{
    this->setX(b.getX());
    this->setY(b.getY());
    this->setZ(b.getZ());
}

R3 R3::operator+(const R3& b)
{
    R3 res(x+b.getX(),y+b.getY(),z+b.getZ());
    return res;
}
R3 R3::operator-(const R3& b)
{
    R3 res(x-b.getX(),y-b.getY(),z-b.getZ());
    return res;
}

R3 R3::operator*(const float& c)
{
    R3 res(c*x,c*y,c*z);
    return res;
}

void R3::operator +=(const R3& b)
{
    this->setX(x+b.getX());
    this->setY(y+b.getY());
    this->setZ(y+b.getZ());
}

void R3::operator -=(const R3& b)
{
    this->setX(x-b.getX());
    this->setY(y-b.getY());
    this->setZ(z-b.getZ());
}

void R3::operator *=(const float& c)
{
    this->setX(c*x);
    this->setY(c*y);
    this->setZ(c*z);
}

std::tuple<float,float,float> R3::operator[](const R3&)
{
    return this->getCoords();
}
