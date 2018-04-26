
/*********************************

**********************************/


#include <QtTest/QtTest>
#include "../C3.h"
#include "../R3.h"

class Tests:public QObject
{
    Q_OBJECT
private slots:
    //C3 Unit Tests
    void C3_Add_1();
    void C3_Add_2();
    void C3_Sub_1();
    void C3_Sub_2();
    void C3_Scal_Mul_1();
    void C3_Scal_Mul_2();
    void C3_Mul_1();
    void C3_Scal_Div_1();
    void C3_Scal_Div_2();
    void C3_Div_1();
    void C3_PS_1();
    void C3_PS_2();
    void C3_PS_3();
    void C3_PS_4();
    void C3_Norm_1();
    void C3_Norm_2();
    void C3_V_Norm_1();
    void C3_V_Norm_2();
    //R3 Unit Tests

};

void Tests::C3_Add_1()
{
    C3 a(0,0,0,1,1,1);
    C3 b(1,1,1,0,0,0);
    QCOMPARE(a+b,C3(1.,1.,1.,1.,1.,1.));
}
void Tests::C3_Add_2()
{
    C3 a(C(0,1),C(1,0),C(1,1));
    C3 b(C(1,0),C(1,1),C(1,1));
    QCOMPARE(a+b,C3(1.,2.,2.,1.,1.,2.));
}
void Tests::C3_Sub_1()
{
    C3 a(0,0,0,1,1,1);
    C3 b(1,1,1,0,0,0);
    QCOMPARE(a-b,C3(-1.,-1.,-1.,1.,1.,1.));
}
void Tests::C3_Sub_2()
{
    C3 a(C(0,1),C(1,0),C(1,1));
    C3 b(C(1,0),C(1,1),C(1,1));
    QCOMPARE(a-b,C3(-1.,0.,0.,1.,-1.,0.));
}
void Tests::C3_Scal_Mul_1()
{
    C3 a(C(1,1),C(1,0),C(1,1));
    R b(7);
    QCOMPARE(a*b,C3(7,7,7,7,0,7));
}
void Tests::C3_Scal_Mul_2()
{
    C3 a(C(1,1),C(1,0),C(1,1));
    C b(7,2);
    QCOMPARE(a*b,C3(5,7,5,9,2,9));
}
void Tests::C3_Mul_1()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    C3 b(C(1,-3),C(2,5),C(0,1));
    QCOMPARE(a*b,C3(C(4,-2),C(2,5),C(2,1)));
}
void Tests::C3_Scal_Div_1()
{
    C3 a(C(1,1),C(1,0),C(1,1));
    R b(7);
    QCOMPARE(a/b,C3(1/7.,1/7.,1/7.,1/7.,0,1/7.));
}
void Tests::C3_Scal_Div_2()
{
    C3 a(C(1,1),C(1,0),C(1,1));
    C b(7,2);
    QCOMPARE(a/b,C3(9/53.,7/53.,9/53.,5/53.,-2/53.,5/53.));
}
void Tests::C3_Div_1()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    C3 b(C(1,-3),C(2,5),C(0,1));
    QCOMPARE(a/b,C3(C(-1/5.,2/5.),C(2/29.,-5/29.),C(-2,-1)));
}
void Tests::C3_PS_1()
{
    C3 a(C(1,1),C(1,1),C(1,1));
    C3 b(C(1,1),C(1,1),C(1,1));
    QCOMPARE((a,b),C(6,0));
}
void Tests::C3_PS_2()
{
    C3 a(C(1,1),C(1,1),C(1,1));
    C3 b(C(1,1),C(1,1),C(1,1));
    QCOMPARE((b,a),C(6,0));
}
void Tests::C3_PS_3()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    C3 b(C(1,-3),C(2,5),C(0,1));
    QCOMPARE((a,b),C(-2,-2));
}
void Tests::C3_PS_4()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    C3 b(C(1,-3),C(2,5),C(0,1));
    QCOMPARE((b,a),C(-2,2));
}
void Tests::C3_Norm_1()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    QCOMPARE(a.norm(),sqrt(8));
}
void Tests::C3_Norm_2()
{
    C3 a(C(1,-3),C(2,5),C(0,1));
    QCOMPARE(a.norm(),2*sqrt(10));
}
void Tests::C3_V_Norm_1()
{
    C3 a(C(1,1),C(1,0),C(1,-2));
    QCOMPARE(a.v_norm(),R3(sqrt(2),1,sqrt(5)));
}
void Tests::C3_V_Norm_2()
{
    C3 a(C(1,-3),C(2,5),C(0,1));
    QCOMPARE(a.v_norm(),R3(sqrt(10),sqrt(29),1));
}

QTEST_MAIN(Tests)
#include "testmodule.moc"

