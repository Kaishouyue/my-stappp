/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with Cbar       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, May 25, 2025                                             */
/*****************************************************************************/



#include "Q4.h"
#include "MatrixMultiplication.h"

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

//    Calculate element stiffness matrix
//    Upper triangular matrix, stored as an array column by column starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix, unsigned int Intmode)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    double ksi=0;
    double eta=0;
    double weight =0;        
    int J_1=1;
    myMatrix adJ(2,2);
    myMatrix B(3,8);
    myMatrix Bt(8,3);
    myMatrix D(3,3);
    myMatrix G(2,4); 
    myMatrix numdaN(2,4);
    myMatrix K(8,8);
    int idx = 0; // Index for storing the stiffness matrix
    double x[4], y[4]; // Coordinates of the 4 nodes
    for (unsigned int i = 0; i < 4; i++) {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }
     //calculate D matrix
        CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
        D(0,0) = material_->E / (1 - material_->nu * material_->nu);
        D(0,1) = material_->nu * D(0,0);
        D(1,0) = D(0,1);
        D(1,1) = D(0,0);   
        D(2,2) = 0.5 * material_->E / (1 + material_->nu);

    // Gauss integration
    if (Intmode == 1) { // 1-point Gauss integration
        ksi = GaussPoints1[0][0];
        eta = GaussPoints1[0][1];
        weight = GaussWeights1[0];
        K=K+ GaussInt4Q4(D, ksi, eta, weight, x, y);
    }
    if (Intmode == 0) { // 4-point Gauss integration
        for (int i = 0; i < 4; ++i) {
            ksi = GaussPoints2[i][0];
            eta = GaussPoints2[i][1];
            weight = GaussWeights2[i];
            K = K + GaussInt4Q4(D, ksi, eta, weight, x, y);
        }
    }
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i <= j; ++i) {
                Matrix[idx++] = K(i, j);
            }
        }   
}

