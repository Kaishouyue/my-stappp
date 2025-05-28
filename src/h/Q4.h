/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with Cbar       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, May 25, 2025                                             */
/*****************************************************************************/

#pragma once
#include "Element.h"
#include "gaussIntegral.h"

using namespace std;

//! Q4 element class
class CQ4 : public CElement
{
public:

//!	Constructor
    CQ4();

//!	Destructor
    ~CQ4();

//!	Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
    virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix, unsigned int Intmode);

//!	Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement);
};