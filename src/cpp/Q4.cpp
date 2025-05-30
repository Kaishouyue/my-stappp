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
    myMatrix D(3,3);
    myMatrix K(8,8);
    int idx = 0; // Index for storing the stiffness matrix
    double x[4], y[4]; // Coordinates of the 4 nodes
    for (unsigned int i = 0; i < 4; i++) {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }
     //calculate D matrix
        CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
        D(0,0) = material_->E / (1 - material_->nu * material_->nu)* material_->t; // E/(1-nu^2) * t
        D(0,1) = material_->nu * D(0,0)* material_->t; // nu * E/(1-nu^2) * t
        D(1,0) = D(0,1);
        D(1,1) = D(0,0);   
        D(2,2) = 0.5 * material_->E / (1 + material_->nu)* material_->t; // 0.5 * E/(1+nu) * t

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
        for (int i = j; i >= 0; --i) {
            Matrix[idx++] = K(i, j);
        }
    }  
}

//    Calculate element stress
void CQ4::ElementStress(double* stress, double* Displacement)
{
    // Calculate the displacement vector for the 4 nodes
   double x[4], y[4]; // Coordinates of the 4 nodes
    for (unsigned int i = 0; i < 4; i++) {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }
    myMatrix u(8,1);
    for (unsigned int i = 0; i < 4; i++) {
        u(2*i, 0) = Displacement[nodes_[i]->bcode[0]];
        u(2*i+1, 0) = Displacement[nodes_[i]->bcode[1]];
    }
    myMatrix D(3,3);
    CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
    D(0,0) = material_->E / (1 - material_->nu * material_->nu)* material_->t; // E/(1-nu^2) * t
    D(0,1) = material_->nu * D(0,0)* material_->t; // nu * E/(1-nu^2) * t
    D(1,0) = D(0,1);
    D(1,1) = D(0,0);   
    D(2,2) = 0.5 * material_->E / (1 + material_->nu)* material_->t; // 0.5 * E/(1+nu) * t
    myMatrix adJ(2,2);
    myMatrix B(3,8);
    myMatrix G(2,4);
    myMatrix numdaN(2,4);
    double J_1 = 1;
    //calculate Jacobian matrix
    adJ(1,1) = 0.25 * ((x[1] - x[0])*1 - (x[3] - x[2])*1);
    adJ(1,0) = -0.25 * ((x[3] - x[0])*1 + (x[2] - x[1])*1);
    adJ(0,1) = -0.25 * ((y[1] - y[0])*1 - (y[3] - y[2])*1);
    adJ(0,0) = 0.25 * ((y[3] - y[0])*1 + (y[2] - y[1])*1);
    J_1 = 1/(adJ(0,0) * adJ(1,1) - adJ(0,1) * adJ(1,0));
        G(0,0) = -1 ;
        G(0,1) = 1 ;
        G(0,2) = 1 ;
        G(0,3) = -1 ;
        G(1,0) = -1 ;
        G(1,1) = -1 ;
        G(1,2) = 1 ;
        G(1,3) = 1 ;
        G = G * 0.25;
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
        myMatrix stress_vector = D * (B * u);
    stress[0] = stress_vector(0, 0); // sigma_xx
    stress[1] = stress_vector(1, 0); // sigma_yy
    stress[2] = stress_vector(2, 0); // sigma_xy
}

    void CQ4::GenerateLocationMatrix()
    {
        unsigned int i = 0;
        for (unsigned int N = 0; N < NEN_; N++)
            for (unsigned int D = 0; D < 2; D++)
                LocationMatrix_[i++] = nodes_[N]->bcode[D];
    }