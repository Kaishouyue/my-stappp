/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with Cbar       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, May 25, 2025                                             */
/*****************************************************************************/



#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//    Constructor
CQ4::CQ4()
{
    NEN_ = 4;    // Each element has 4 nodes
    nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];
    
    ElementMaterial_ = nullptr;
}

//    Destructor
CQ4::~CQ4()
{
}

//    Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;    // Material property set number
    unsigned int N1, N2, N3, N4;    // Node number of 4 nodes
    
    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
    
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    
    return true;
}

//    Write element data to stream
void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}