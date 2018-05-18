#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include "RNM/RNM.hpp"
#include "C3.h"
#include "maillage3d.h"
//Sparse matrix using the Morse decompsition by compressed row

class VectorC3
{
public:
    KN<C> X[3];
    VectorC3(){}

    ~VectorC3(){X[0].destroy();X[1].destroy();X[2].destroy();}

    VectorC3(KN<C> x,KN<C> y,KN<C> z){X[0]=x;X[1]=y;X[2]=z;}
    VectorC3(ifstream& file);

    int size(){
        if( X[0].size()==X[1].size() && X[1].size()==X[2].size())
            return X[0].size();
        else
            return -1;}
    void resize(long nn){X[0].resize(nn);X[1].resize(nn);X[2].resize(nn);}
    void init(long nn){X[0].init(nn);X[1].init(nn);X[2].init(nn);}
    R Re_norm(int _i=-1)
    //If an argument ins imputed, the method will only compute the
    //norm of the imputed component, else it will compute the entire norm of the vector
    {
        //Definition of the norm of a (complex vectorial) vector i.e. #C^(dxn)
        //                 _                                 _
        //                |   3    N                          |
        // || V || = Sqrt |  Sum [Sum <V.X[i][j],V.X[i][j]> ] |
        //                |_ i=1  j=1                        _|
        int i;
        R sq_res=0.;
        if (_i==-1)
        {
            for (i=0;i<3;i++)
            {
                sq_res+= (this->X[i].norm());
            }
            return sqrt(sq_res);
        }
        else
            return this->X[_i].l2();
    }

    C3 Re_norms()
    //Pour ne pas avoir de souci avec les variables issues de produits scalaires
    //ou de normes, on suppose que la norme d'un vecteur C3 est une donnée de C3
    // meme si sa partie imaginaire est évidemment nulle.
    {return C3(X[0].l2(),X[1].l2(),X[2].l2(),0,0,0);}

    //Operateurs d'affectation (+,-,*,/)
    VectorC3& operator+=(VectorC3 P) {
        X[0] += P.X[0];
        X[1] += P.X[1];
        X[2] += P.X[2];
        return  *this;}

    VectorC3& operator-=(VectorC3 P) {
        X[0] -= P.X[0];
        X[1] -= P.X[1];
        X[2] -= P.X[2];
        return  *this;}

    VectorC3& operator*=(const R k) {
        X[0] *= k;
        X[1] *= k;
        X[2] *= k;
        return  *this;}

    VectorC3& operator*=(const C k) {
        X[0] *= k;
        X[1] *= k;
        X[2] *= k;
        return  *this;}

    VectorC3& operator*=(C3 k) {
        X[0] *= k.X_();
        X[1] *= k.Y_();
        X[2] *= k.Z_();
        return  *this;}

    VectorC3& operator*=(R3 k) {
        X[0] *= k.X_();
        X[1] *= k.Y_();
        X[2] *= k.Z_();
        return  *this;}

    VectorC3& operator/=(const R k) {
        X[0] /= k;
        X[1] /= k;
        X[2] /= k;
        return  *this;}

    VectorC3& operator/=(R3 k) {
        X[0] /= k.X_();
        X[1] /= k.Y_();
        X[2] /= k.Z_();
        return  *this;}

    VectorC3& operator/=(C3 k) {
        X[0] /= k.X_();
        X[1] /= k.Y_();
        X[2] /= k.Z_();
        return  *this;}

    VectorC3& operator/=(const C k) {
        X[0] /= k;
        X[1] /= k;
        X[2] /= k;
        return  *this;}

    void operator=(VectorC3 P) {
        this->X[0].resize(P.size()); this->X[0]=P.X[0];
        this->X[1].resize(P.size()); this->X[1]=P.X[1];
        this->X[2].resize(P.size()); this->X[2]=P.X[2];
        }

    VectorC3 operator+(VectorC3 P) {
        KN<C> x(X[0]);x.resize(P.size()); x+=P.X[0];
        KN<C> y(X[1]);y.resize(P.size()); y+=P.X[1];
        KN<C> z(X[2]);z.resize(P.size()); z+=P.X[2];
        return VectorC3(x,y,z);}

    VectorC3 operator-(VectorC3 P) {
        KN<C> x(X[0]);  x.resize(P.size());x-=P.X[0];
        KN<C> y(X[1]);  y.resize(P.size());y-=P.X[1];
        KN<C> z(X[2]);  z.resize(P.size());z-=P.X[2];
        return VectorC3(x,y,z);}

    VectorC3 operator*(C k) {
        KN<C> x(X[0]); x.resize(this->size()); x*=k;
        KN<C> y(X[1]); y.resize(this->size());y*=k;
        KN<C> z(X[2]); z.resize(this->size());z*=k;
        return VectorC3(x,y,z);}

    VectorC3 operator*(const R3 & k) {
        KN<C> x(X[0]); x.resize(this->size()); x*=k.X_();
        KN<C> y(X[1]); y.resize(this->size());y*=k.Y_();
        KN<C> z(X[2]); z.resize(this->size());z*=k.Z_();
        return VectorC3(x,y,z);}

    VectorC3 operator*(const C3 & k) {
        KN<C> x(X[0]); x.resize(this->size()); x*=k.X_();
        KN<C> y(X[1]); y.resize(this->size());y*=k.Y_();
        KN<C> z(X[2]); z.resize(this->size());z*=k.Z_();
        return VectorC3(x,y,z);}

    C3 operator, (VectorC3 P)
    {
        C dp_0,dp_1,dp_2;
        dp_0 = (KN<C>(X[0]),(conj_KN_<C>(P.X[0])));
        dp_1 = (KN<C>(X[1]),(KN<C>(P.X[1])));
        dp_2 = (KN<C>(X[2]),(KN<C>(P.X[2])));
        return(C3(dp_0,dp_1,dp_2));}
};

typedef VectorC3 MatDiag;

class MatSparseC3
{
public:
    MatSparseC3();
    MatSparseC3(ifstream& file);

    int DiagMatMul(MatDiag);
    int MatDiagMul(MatDiag);

    void MatVectMul (VectorC3 &i,VectorC3 &o,int _comp=-1)const;
    int EditProfileFromMesh(Maillage3D);
    int Insert(int,int);

    int IdxSearch( int i,int j);

    VectorC3 &GetData(){return C_data;}
    KN<int> &GetJ() {return J;}
    KN<int> &GetP() {return P;}
    int GetSize() {return P.size();}
private:
    VectorC3 C_data;
    KN<int> J;
    KN<int> P;
    int size;
};

#endif // SPARSEMAT_H
