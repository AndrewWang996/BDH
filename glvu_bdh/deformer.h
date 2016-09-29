#pragma once

#include "vaomesh.h"

std::string datadir;
int numMeshVertex = 20000;

struct Deformer
{
    bool needIteration = false;
    virtual std::string name(){ return "UNKNOWN"; }
    virtual void preprocess() {}
    virtual void updateP2PConstraints(int){}
    virtual void deform() = 0;
    virtual void resetDeform(){}
    virtual void getResult(){}
    virtual void saveData(){}
};

