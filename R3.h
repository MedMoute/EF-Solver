#ifndef R3_H
#define R3_H

#include <cassert>

using namespace std;
//Toutes les fonctions sont implémentées en inline dans ce header.

typedef double R; //Assimillation des réels aux double

class R3
{
public:
    //Constructeur par initialisation à zero
    R3():X(0.),Y(0.),Z(0.){}
    //Constructeur explictite à partir de 3 réels
    R3(R a,R b,R c):X(a),Y(b),Z(c){}

    //Operateurs d'affectation (+,-,*,/)
    R3& operator+=(const R3 & P) {
        X += P.X;
        Y += P.Y;
        Z += P.Z;
        return  *this;}

    R3& operator-=(const R3 & P) {
        X -= P.X;
        Y -= P.Y;
        Z -= P.Z;
        return  *this;}

    R3& operator*=(const R & k) {
        X *= k;
        Y *= k;
        Z *= k;
        return  *this;}

    R3& operator/=(const R & k) {
        X /= k;
        Y /= k;
        Z /= k;
        return  *this;}


    bool operator==(R3 P)const{
        if (X==P.X_() && Y==P.Y_() && Z==P.Z_())
            return true;
        else return false;
    }

    //Opérateurs binaires (+,-,*,/,PS(,))
    R3 operator+(const R3 & P) {
        return(R3(X+P.X_(),Y+P.Y_(),Z+P.Z_()));}
    R3 operator-(const R3 & P) {
        return(R3(X-P.X_(),Y-P.Y_(),Z-P.Z_()));}
    R3 operator*(const R & k) {
        return(R3(k*X,k*Y,k*Z));}
    R3 operator*(const R3 & P){
        return(R3(X*P.X_(),Y*P.Y_(),Z*P.Z_()));}
    R3 operator/(const R & k){
        return(R3(k/X,k/Y,k/Z));}
    R3 operator/(const R3 & P){
        return(R3(X/P.X_(),Y/P.Y_(),Z/P.Z_()));}
    R operator, (const R3 & P) {
        return(R(X*P.X_()+Y*P.Y_()+Z*P.Z_()));}

    //Opérateurs unaires(+(Identité),-(Négatif),*(Conjugate))
    R3 operator+()const {return *this;}

    R3 operator-()const {return R3(-X,-Y,-Z);}

    R3 sqrt_() const{return R3(sqrt(X),sqrt(Y),sqrt(Z));}

    //Operateur tableau : récupère les composantes dans l'ordre de construction
    R & operator[](int i) { //Operateur autorisant les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else if (i==2) return Z;
        else
        {
            assert(0);
            exit(1);
        }
    }

    const R & operator[](int i) const{ //Operateur n'autorisant pas les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else if (i==2) return Z;
        else
        {
            assert(0);
            exit(1);
        }
    }

    //Méthodes

    R3 div_(const R3& P)const {
        return(R3(X/P.X_(),Y/P.Y_(),Z/P.Z_()));
    }
    //Récupération des composantes
    R comp(int i) const{
        if (i==1) return (R(X));
        else if (i==2) return (R(Y));
        else if (i==3) return (R(Z));
        else
        {
            assert(0);
            exit(1);
        }
    }
    //Méthode pour obtenir les coefficients inverses
    R3 inv() const{return R3(1/X,1/Y,1/Z);}

    double norm() const{return sqrt(X*X+Y*Y+Z*Z);}

    R X_()const{return X;}
    R Y_()const{return Y;}
    R Z_()const{return Z;}
private:
    //Membres
    R X,Y,Z; //on interdit l'accès direct aux données en dehors de la classe.
};
inline ostream & operator<<(ostream & f, R3 c)
{
f << c.X_() <<" // "<< c.Y_()<<" // "<< c.Z_() ;
return f;
}
#endif // R3_H
