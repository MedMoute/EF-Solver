#ifndef UTILS_H
#define UTILS_H

#include<cmath>
#include<algorithm>
#include<tuple>

#include "R2.h"
#include "R3.h"

using namespace std;

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
};

#endif // UTILS_H
