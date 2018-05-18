#include <iostream>

#include "utils.h"


float util::euclDistance( R2& a , R2& b)
{
    float dx = max(a.X_()-b.X_(),b.X_()-a.X_());
    float dy = max(a.Y_()-b.Y_(),b.Y_()-a.Y_());
    return (sqrt(dx*dx + dy*dy));
}

float util::euclDistance( R3& a , R3& b)
{
    float dx = max(a.X_()-b.X_(),b.X_()-a.X_());
    float dy = max(a.Y_()-b.Y_(),b.Y_()-a.Y_());
    float dz = max(a.Z_()-b.Z_(),b.Z_()-a.Z_());
    return (sqrt(dx*dx + dy*dy + dz*dz));
}

float util::scalPdt( R2& a,  R2& b)
{
    return (a.X_()*b.X_() + a.Y_()*b.Y_());
}

float util::scalPdt( R3& a,  R3& b)
{
     return (a.X_()*b.X_() + a.Y_()*b.Y_()+a.Z_()*b.Z_());
}

// correspondance entre le type d'element et le nb de sommets de cet element
int util::type2nnodes(int type) {

    int nnodes[31] = {2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56} ;
    return nnodes[type-1] ;

}

void util::print_separator()
{
    cout<<"-----------------------"<<endl;
}

int util::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int util::kron(int i,int j)
{
    return (i==j)?1:0;
}
double util::simplex3Mesure(R3 a,R3 b,R3 c,R3 d)
{
    R3 p0_1=b-a;
    R3 p0_2=c-a;
    R3 p0_3=d-a;

    double m1=p0_1.X_()*p0_2.Y_()*p0_3.Z_()-(p0_1.Z_()*p0_2.Y_()*p0_3.X_());
    double m2=p0_1.Y_()*p0_2.Z_()*p0_3.X_()-(p0_1.X_()*p0_2.Z_()*p0_3.Y_());
    double m3=p0_1.Z_()*p0_2.X_()*p0_3.Y_()-(p0_1.Y_()*p0_2.X_()*p0_3.Z_());
    return ((m1+m2+m3)/6.);
}

double util::simplex2Mesure(R3 a,R3 b,R3 c)
{
    //Using Heron formula
    R3 p0_1=b-a;
    R3 p0_2=c-a;
    R3 p1_2=c-b;
    double d0_1=p0_1.norm();
    double d0_2=p0_2.norm();
    double d1_2=p1_2.norm();

    double s=(d0_1+d0_2+d1_2)/2;

    return(sqrt(s*(s-d0_1)*(s-d0_2)*(s-d1_2)));
}
