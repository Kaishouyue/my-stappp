/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E>>  nu  >>  t;	// Young's modulus , Poisson's ratio and thickness
    if (nu < 0.0 || nu > 0.5)
	{
		cerr << "*** Error *** Poisson's ratio must be between 0 and 0.5 !" << endl
			 << "    Provided value : " << nu << endl;
		return false;
	}
	return true;
}

//	Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << t << endl;
}
//	Material class for 4Q element