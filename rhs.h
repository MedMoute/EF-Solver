//Class implementing the Right Hand Side data parsing for
//the MAP program
#ifndef RHS_H
#define RHS_H

#include<vector>
#include <fstream>
#include <iostream>
#include <map>

#include "C3.h"
#include "muparserx/parser/mpParser.h"

using namespace mup;

class RHS
{
public:
    //Reading the RHS Data from a file
    RHS(ifstream& file, map<int,R3>* node_data);

    string* GetExpr(){return expr;}
    vector<C3>* GetRHS(){return rhs;}
    int GetSize(){return rhs->size();}
private:
    vector<C3>* rhs;
    string expr[3];

};

#endif // RHS_H
