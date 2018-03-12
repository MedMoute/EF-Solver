#include "r2.h"

R2::R2()
{
    x=0;
    y=0;
}

R2::R2(float _x, float _y)
{
    x=_x;
    y=_y;
}

R2::R2(const R2& point)
{
    std::pair<float,float> data=point.getCoords();
    x=data.first;
    y=data.second;
}

std::pair<float,float> R2::getCoords() const
{
    return std::make_pair(x,y);
}

float R2::getX() const
{
    return x;
}

float R2::getY() const
{
    return y;
}

void R2::setX(float _x)
{
    x=_x;
}

void R2::setY(float _y)
{
    y=_y;
}

void R2::displayCoords()
{
    std::cout<< "Point data : x ="<<getX()<<" --- y = "<<getY()<<std::endl;
}

bool R2::operator ==(const R2& b)
{
    return ((b.getX()==x)&&(b.getY()==y));
}

std::ostream& R2::operator<<(std::ostream& out)
{
    out << "(" << x << ", " << y << ")";
    return out;
}

void R2::operator =(const R2& b)
{
    this->setX(b.getX());
    this->setY(b.getY());
}

R2 R2::operator+(const R2& b)
{
    R2 res(x+b.getX(),y+b.getY());
    return res;
}
R2 R2::operator-(const R2& b)
{
    R2 res(x-b.getX(),y-b.getY());
    return res;
}

R2 R2::operator*(const float& c)
{
    R2 res(c*x,c*y);
    return res;
}

void R2::operator +=(const R2& b)
{
    this->setX(x+b.getX());
    this->setY(y+b.getY());
}

void R2::operator -=(const R2& b)
{
    this->setX(x-b.getX());
    this->setY(y-b.getY());
}

void R2::operator *=(const float& c)
{
    this->setX(c*x);
    this->setY(c*y);
}

std::pair<float,float> R2::operator[](const R2&)
{
    return this->getCoords();
}
