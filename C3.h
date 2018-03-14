#ifndef C3_H
#define C3_H

#include <complex>
#include <cassert>

#include "R3.h"

using namespace std;
//Définition de la classe C3 : Nombres complexes en dimension 3
//Toutes les fonctions sont implémentées en inline dans ce header.

typedef double R; //Assimillation des réels aux double
typedef complex<double> C;


class C3
{
public:

    //Constructeur par initialisation à zero
    C3():X(0.,0.),Y(0.,0.),Z(0.,0.){}
    //Constructeur explictite à partir de 6 réels
    C3(R a,R b,R c,R d,R e, R f):X(a,d),Y(b,e),Z(c,f){}
    //Constructeur explicite à partir de 3 complexes
    C3(C x,C y,C z):X(x.real(),x.imag()),Y(y.real(),y.imag()),Z(z.real(),z.imag()){}
    //Constructeur explicite à partir de 2 elements de R3
    C3(R3 a,R3 b):X(a.X_(),b.X_()),Y(a.Y_(),b.Y_()),Z(a.Z_(),b.Z_()){}

    //Operateurs d'affectation (+,-,*,/)
    C3& operator+=(C3 P) {
        X += P.X_();
        Y += P.Y_();
        Z += P.Z_();
        return  *this;}

    C3& operator-=(C3 P) {
        X -= P.X_();
        Y -= P.Y_();
        Z -= P.Z_();
        return  *this;}

    C3& operator*=(const R k) {
        X *= k;
        Y *= k;
        Z *= k;
        return  *this;}

    C3& operator/=(const R k) {
        X /= k;
        Y /= k;
        Z /= k;
        return  *this;}


    C3& operator*=(C3 P) {
        X *= P.X_();
        Y *= P.Y_();
        Z *= P.Z_();
        return  *this;}

    C3& operator/=(C3 P) {
        X /= P.X_();
        Y /= P.Y_();
        Z /= P.Z_();
        return  *this;}

    //Opérateurs binaires (+,-,*,/,PS(,))
    C3 operator+(C3 P) {
        return(C3(X.real()+P.X_().real(),Y.real()+P.Y_().real(),Z.real()+P.Z_().real(),X.imag()+P.X_().imag(),Y.imag()+P.Y_().imag(),Z.imag()+P.Z_().imag()));}
    C3 operator-(C3 P) {
        return(C3(X.real()-P.X_().real(),Y.real()-P.Y_().real(),Z.real()-P.Z_().real(),X.imag()-P.X_().imag(),Y.imag()-P.Y_().imag(),Z.imag()-P.Z_().imag()));}
    C3 operator*( R k) {
        return(C3(k*X.real(),k*Y.real(),k*Z.real(),k*X.imag(),k*Y.imag(),k*Z.imag()));}
    C3 operator* (C3 P) {
        return(C3(C(X.real()*P.X_().real()-X.imag()*P.X_().imag(),X.imag()*P.X_().real()+X.real()*P.X_().imag()),C(Y.real()*P.Y_().real()-Y.imag()*P.Y_().imag(),Y.imag()*P.Y_().real()+Y.real()*P.Y_().imag()),C(Z.real()*P.Z_().real()-Z.imag()*P.Z_().imag(),Z.imag()*P.Z_().real()+Z.real()*P.Z_().imag())));}

    C3 operator/(const R k){
        return(C3((1/k)*X.real(),(1/k)*Y.real(),(1/k)*Z.real(),(1/k)*X.imag(),(1/k)*Y.imag(),(1/k)*Z.imag()));}
  C3 operator/ (C3 P){
        return(C3(X.real()/P.v_norm().X_(),Y.real()/P.v_norm().Y_(),Z.real()/P.v_norm().Z_(),X.imag()/P.v_norm().X_(),Y.imag()/P.v_norm().Y_(),Z.imag()/P.v_norm().Z_()));}

    // Attention le produit scalaire dans C3 (donc hermitien) n'est PAS commutatif
    //          3        ___      3
    // <x,y> = Sum [ x_i*y_i ] = Sum [Re(x_i)*Re(y_i)+Im(x_i)*Im(y_i)+i*(Im(x_i)*Re(y_i)-Re(x_i)*Im(y_i))]
    //         i=1               i=1
    C operator, (C3 & P) {
        return(C(X.real()*P.X_().real()+X.imag()*P.X_().imag()+Y.real()*P.Y_().real()+Y.imag()*P.Y_().imag()+Z.real()*P.Z_().real()+Z.imag()*P.Z_().imag(),X.imag()*P.X_().real()-X.real()*P.X_().imag()+Y.imag()*P.Y_().real()-Y.real()*P.Y_().imag()+Z.imag()*P.Z_().real()-Z.real()*P.Z_().imag()));}

    //Opérateurs unaires(+(Identité),-(Négatif),*(Conjugate))
    C3 operator+()const {return *this;}

    C3 operator-()const {return C3(-X.real(),-Y.real(),-Z.real(),-X.imag(),-Y.imag(),-Z.imag());}

    C3 operator*()const {return C3(X.real(),Y.real(),Z.real(),-X.imag(),-Y.imag(),-Z.imag());}

    //Operateur tableau : récupère les composantes dans l'ordre de construction
    C & operator[](int i) { //Operateur autorisant les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else if (i==2) return Z;
        else
        {
            assert(0);
            exit(1);
        }
    }

    const C & operator[](int i) const{ //Operateur n'autorisant pas les modification
        if(i==0) return X;
        else if (i==1) return Y;
        else if (i==2) return Z;
        else
        {
            assert(0);
            exit(1);
        }
    }
    R3 v_norm()
    {return R3(sqrt(X.real()*X.real()+X.imag()*X.imag()),sqrt(Y.real()*Y.real()+Y.imag()*Y.imag()),sqrt(Z.real()*Z.real()+Z.imag()*Z.imag()));}
    //Méthodes
    R norm()
    {return (sqrt(X.real()*X.real()+Y.real()*Y.real()+Z.real()*Z.real()+X.imag()*X.imag()+Y.imag()*Y.imag()+Z.imag()*Z.imag()));}

    //Récupération des composantes
    C comp(int i) const{
        if (i==1) return (C(X.real(),X.imag()));
        else if (i==2) return (C(Y.real(),Y.imag()));
        else if (i==3) return (C(Z.real(),Z.imag()));
        else
        {
            assert(0);
            exit(1);
        }
    }
    //Recuperation de la partie reele
    R3 Re() const{ return (R3(X.real(),Y.real(),Z.real())); }
    //Recupération de la partie imaginaire
    R3 Im() const{ return (R3(X.imag(),Y.imag(),Z.imag())); }

    C X_(){return X;}
    C Y_(){return Y;}
    C Z_(){return Z;}

private:
    //Membres
    complex<R> X,Y,Z; //on interdit l'accès direct aux données en dehors de la classe.
};

#endif // C3_H
