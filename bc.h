//Class implementing the Boundary Conditions data parsing for
//the MAP program

#ifndef BC_H
#define BC_H

#include<vector>
#include <fstream>
#include <iostream>
#include <map>

#include "C3.h"

class BC
{
public:
    //Reading the RHS Data from a file
    BC(ifstream& file, map<int,R3>* node_data);

    int GetSize(){return bc->size();}
    string* GetExpr(){return expr;}
    vector<C3>* GetBC(){return bc;}

private:
    vector<C3>* bc;
    string expr[3];
};

#endif // BC_H
