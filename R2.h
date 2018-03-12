#ifndef R2_H
#define R2_H

#include <cassert>
using namespace std;
//Toutes les fonctions sont implémentées en inline dans ce header.

typedef double R; //Assimillation des réels aux double

class R2
{
public:
    //Constructeur par initialisation à zero
    R2():X(0.),Y(0.){}
    //Constructeur explictite à partir de 2 réels
    R2(R a,R b):X(a),Y(b){}

    //Operateurs d'affectation (+,-,*,/)
    R2& operator+=(const R2 & P) {
        X += P.X;
        Y += P.Y;
        return  *this;}

    R2& operator-=(const R2 & P) {
        X -= P.X;
        Y -= P.Y;
        return  *this;}

    R2& operator*=(const R & k) {
        X *= k;
        Y *= k;
        return  *this;}

    R2& operator/=(const R & k) {
        X /= k;
        Y /= k;
        return  *this;}

    //Opérateurs binaires (+,-,*,/,PS(,))
    R2 operator+(R2 & P) {
        return(R2(X+P.X_(),Y+P.Y_()));}
    R2 operator-(R2 & P) {
        return(R2(X-P.X_(),Y-P.Y_()));}
    R2 operator*(R & k) {
        return(R2(k*X,k*Y));}
    R2 operator/(R & k){
        return(R2(k/X,k/Y));}

    R operator, (R2 & P) {
        return(R(X*P.X_()+Y*P.Y_()));}

    //Opérateurs unaires(+(Identité),-(Négatif))
    R2 operator+()const {return *this;}

    R2 operator-()const {return R2(-X,-Y);}

    //Operateur tableau : récupère les composantes dans l'ordre de construction
    R & operator[](int i) { //Operateur autorisant les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else
        {
            assert(0);
            exit(1);
        }
    }

    const R & operator[](int i) const{ //Operateur n'autorisant pas les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else
        {
            assert(0);
            exit(1);
        }
    }

    //Méthodes

    //Récupération des composantes
    R comp(int i) const{
        if (i==1) return (R(X));
        else if (i==2) return (R(Y));
        else
        {
            assert(0);
            exit(1);
        }
    }

    R X_(){return X;}
    R Y_(){return Y;}
private:
    //Membres
    R X,Y; //on interdit l'accès direct aux données en dehors de la classe.
};

#endif // R2_H
