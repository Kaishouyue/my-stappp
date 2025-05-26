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
    int ksi=0;
    int eta=0;
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

    if (Intmode == 1){// reduce integral
    //calculate Jacobian matrix
        adJ(1,1) = 0.25 * (x[1] - x[0] - (x[3] - x[2]));
        adJ(1,0) = -0.25 * (x[3] - x[0] + (x[2] - x[1]));
        adJ(0,1) = -0.25 * (y[1] - y[0] - (y[3] - y[2]));
        adJ(0,0) = 0.25 * (y[3] - y[0] + (y[2] - y[1]));
        J_1 = 1/(adJ(0,0) * adJ(1,1) - adJ(0,1) * adJ(1,0));
        G(0,0) = -1 ;
        G(0,1) = 1 ;
        G(0,2) = 1 ;
        G(0,3) = -1 ;
        G(1,0) = -1 ;
        G(1,1) = -1 ;
        G(1,2) = 1 ;
        G(1,3) = 1 ;
        numdaN=(adJ * G) * J_1; // N = J^-1 * G
        //calculate B matrix
        B(0,0) = numdaN(0,0);
        B(0,2) = numdaN(0,1);
        B(0,4) = numdaN(0,2);
        B(0,6) = numdaN(0,3);
        B(1,1) = numdaN(1,0);
        B(1,3) = numdaN(1,1);
        B(1,5) = numdaN(1,2);
        B(1,7) = numdaN(1,3);
        B(2,0) = numdaN(1,0);
        B(2,1) = numdaN(0,0);
        B(2,2) = numdaN(1,1);
        B(2,3) = numdaN(0,1);
        B(2,4) = numdaN(1,2);
        B(2,5) = numdaN(0,2);
        B(2,6) = numdaN(1,3);
        B(2,7) = numdaN(0,3);
        Bt = B.transpose();
        //calculate D matrix
        CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
        D(0,0) = material_->E / (1 - material_->nu * material_->nu);
        D(0,1) = material_->nu * D(0,0);
        D(1,0) = D(0,1);
        D(1,1) = D(0,0);   
        D(2,2) = 0.5 * material_->E / (1 + material_->nu);
        //calculate K matrix
        K = (GAUSS_1_WEIGHT[0]* GAUSS_1_WEIGHT[0])* (Bt * D * B) * J_1; // K = sum  B^T * D * B * J^-1
        //store K matrix to Matrix
        idx = 0;
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i <= j; ++i) {
                Matrix[idx++] = K(i, j);
            }
        }
    }
}

