/*//////////////////////////////////////////////////////////////////////////////////////////////////
///  ELC.h   Clique Reduction by Excludable Local Configuration
///  Version 1.04         September 12th, 2014
////////////////////////////////////////////////////////////////////////////////////////////////////

Copyright 2014 Hiroshi Ishikawa. All rights reserved.
This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Hiroshi Ishikawa).

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For the latest version, check: http://www.f.waseda.jp/hfs/indexE.html

////////////////////////////////////////////////////////////////////////////////////////////////////

Software for minimizing a higher-order function of binary variables x_1,...,x_n.
What it actually does is to reduce the function into a first-order MRF, or a
Quadratic Pseudo-Boolean function, i.e., a function of the form
E(x_1, ..., x_n, ..., x_m) = \sum_i Ei(x_i) + \sum_{ij} Eij(x_i,x_j).
possibly with additional variables.
Once the function is reduced, minimization itself can be done with known first-order methods,
such as the QPBO (roof dual) algorithm below.

In this software, there are two ways to reduce a higher-order monomial to first order (second degree).


A) Finding Excludable Local Configuration (ELC)
Reduces the monomial without adding any variable. For some monomials, it cannot be done.
So the B) below must be used after A). The technique is described in the following paper:

[1] Hiroshi Ishikawa. "Higher-Order Clique Reduction without Auxiliary Variables."
In CVPR2014, Columbus, Ohio. June 23-28, 2014.

Additionally, this software provides an approximation that always reduces the term
without additional variables. It does so correctly most of the time, but there is no guarantee.
It is much faster, though. This is also described in [1].


B) Higher-Order Clique Reduction (HOCR)
Additional variables are added to reduce the order of the energy.
The number of variables increases exponentially as the order of the given energy increases.
The technique is described in the following papers:

[2] Hiroshi Ishikawa. "Transformation of General Binary MRF Minimization to the First Order Case."
IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 33, no. 6, pp. 1234-1249,
June 2011.

[3] Hiroshi Ishikawa. "Higher-Order Clique Reduction in Binary Graph Cut."
In CVPR2009: IEEE Computer Society Conference on Computer Vision and Pattern Recognition,
Miami Beach, Florida. June 20-25, 2009.

////////////////////////////////////////////////////////////////////////////////////////////////////

This software is implemented so that it can be used most conveniently
in combination with the QPBO software by Vladimir Kolmogorov available
at http://pub.ist.ac.at/~vnk/software.html

This software has been tested on Windows 7 (x64) with Visual Studio 2010,
Ubuntu 12.04 with g++ 4.6.3, and Ubuntu 12.10 with g++ 4.8.1.
Any report on bugs and results of trying on other platforms is appreciated.

////////////////////////////////////////////////////////////////////////////////////////////////////

Example usage:
Minimize E(x, y, z, w) = x + 4y - z - 2w(y-1) + xy(z+1) - xw(y+1)(z+1), where x,y,z,w are in {0,1}.
This is a third order (fourth degree) function.

Step 1. Use the Pseudo-Boolean function object to build the energy function object.
Step 2. Convert to quadratic, then convert it to QPBO object.
Step 3. Minimize the QPBO object using the QPBO software.
Step 4. Read out the results.

//-------- Sample code -----------
#include "ELC/ELC.h"
#include "QPBO/QPBO.h"
using namespace ELCReduce;
typedef int REAL;


	// mode 0: reduce only the terms with ELCs, convert the rest with HOCR
	// mode 1: reduce all higher-order terms using the approximation
	// mode 2: use only HOCR
void reduce(const PBF<REAL>& pbf, PBF<REAL>& qpbf, int mode, int newvar)
{
	if (mode == 0)
	{
		PBF<REAL> pbf2 = pbf;
		pbf2.reduceHigher(); // Use the ELC technique to reduce higher-order terms without auxiliary variables
		pbf2.toQuadratic(qpbf, newvar); // Reduce the remaining higher-order terms using HOCR adding auxiliary variables
	}
	else if (mode == 1)
	{
		qpbf = pbf;
		qpbf.reduceHigherApprox(); // Use the approximate ELC technique to reduce higher-order terms without auxiliary variables
	}
	else if (mode == 2)
		pbf.toQuadratic(qpbf, newvar); // Reduce to Quadratic pseudo-Boolean function using HOCR.
}

void main()
{
//  Step 1. Use the Pseudo-Boolean function object to build the energy function object.
	PBF<REAL> pbf;
	pbf.AddUnaryTerm(0, 0, 1); // Add the term x
	pbf.AddUnaryTerm(1, 0, 4); // Add the term 4y
	pbf.AddUnaryTerm(2, 0, -1); // Add the term -z
	pbf.AddPairwiseTerm(1, 3, 0, 2, 0, 0); // Add the term -2(y-1)w
	int vars3[3] = {0,1, 2};
	REAL vals3[8] = {0, 0, 0, 0, 0, 0, 1, 2};
	pbf.AddHigherTerm(3, vars3, vals3); // Add the term  xy(z+1)
	int vars4[4] = {0,1, 2, 3};
	REAL vals4[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -2, 0, -2, 0, -4};
	pbf.AddHigherTerm(4, vars4, vals4); // Add the term  -xw(y+1)(z+1)

//  Step 2. Convert to quadratic, then convert it to QPBO object.
	PBF<REAL> qpbf; // quadratic pbf
	int mode = 0;
	reduce(pbf, qpbf, mode, 4); // see above
	int numvars = qpbf.maxID(); // Number of variables
	QPBO<int> qpbo(numvars, numvars * 4);
	qpbf.convert(qpbo, 4); // copy to QPBO object by V. Kolmogorov.
	qpbf.clear(); // free memory

//	Step 3. Minimize the QPBO object using the QPBO software.
	qpbo.MergeParallelEdges();
	qpbo.Solve();
	qpbo.ComputeWeakPersistencies();

//	Step 4. Read out the results.
	int x = qpbo.GetLabel(0);  // negative number means "not sure"
	int y = qpbo.GetLabel(1);
	int z = qpbo.GetLabel(2);
	int w = qpbo.GetLabel(3);
	printf("Solution: x=%d, y=%d, z=%d, w=%d\n", x, y, z, w);
}
//-------- End sample code --------

//////////////////////////////////////////////////////////////////////////////////////////////////*/


#ifndef __ELC_H__
#define __ELC_H__

#include<vector>
#include<algorithm>
#include <iterator>
#include<assert.h>
#include <limits>
#include <iostream>
#include <cstdio>

namespace ELCReduce {

typedef int VID;	// Type for the variable index
typedef std::vector<VID> VVec;
typedef VVec::iterator VVecIt;
typedef unsigned int BITMASK;
const int MASKSIZE = 8 * sizeof(BITMASK);
const int MAXD = MASKSIZE; // Maximum degree allowed
const int MAXVARELC = 16; // Maximum number of variables in exhaustive ELC search

#include "ELC0.h" // internal classes and functions


////////////////////////////////////////////////////////////////////////////////////////////////////
// class PBF: Pseudo-Boolean Function
// Represents a pseudo-Boolean function. Includes a reduction to quadratic pbf.
// REAL: type for coefficients

template <typename REAL>
class PBF
{
	typedef typename Terms<REAL>::TermPointer TermPointer;
	typedef typename Terms<REAL>::variables variables;
public:
	PBF(int maxVars = 0) : constant(0), terms(maxVars) {}

	// clears to free memory
	void clear() {terms.clear();}

	//Print terms
	void print() {terms.print(constant);}

    void printStar()
	{
	    int j = 0;
		for (iterator it = begin(); it != end(); ++it)
		{
			int d = it.degree();
			if (d < 3)
				break;
			smallPBF<REAL> pbf;
			getStar(pbf, it);
			std::printf("Star:\n");
			auto idx = it.vars();
            std::printf("Variables: %i, %i, %i, %i\n", *idx, *(idx + 1), *(idx + 2), *(idx + 3));
			pbf.print();
			if (j==2)
                break;
            j++;
		}
	}

	void printCore()
	{
	    int j = 0;
		for (iterator it = begin(); it != end(); ++it)
		{
			int d = it.degree();
			if (d < 3)
				break;
			smallPBF<REAL> pbf;
			getCore(pbf, it);
            std::printf("Core:\n");
            auto idx = it.vars();
            std::printf("Variables: %i, %i, %i, %i\n", *idx, *(idx + 1), *(idx + 2), *(idx + 3));
			pbf.print();
			if (j==1)
                break;
            j++;
		}
	}

	// Adds unary term Ei(x_i) to the energy function with cost values Ei(0)=E0, Ei(1)=E1.
	// This adds the terms  E0(1 - x_i) + E1 x_i
	void AddUnaryTerm(VID i, REAL E0, REAL E1) {constant += E0; add1(E1 - E0, i);}

	// Adds pairwise term Eij(x_i, x_j) with cost values E00, E01, E10, E11.
	// This adds the terms  E00(1-x_i)(1-x_j) + E01(1-x_i)x_j + E10 x_i(1-x_j) + E11 x_i x_j
	// i < j must hold.
	void AddPairwiseTerm(VID i, VID j, REAL E00, REAL E01, REAL E10, REAL E11)
	{
		constant += E00;
		add1(E10 - E00, i);
		add1(E01 - E00, j);
		add2(E00 - E01 - E10 + E11, i, j);
	}

	// Adds higher-order term Eij...k(x_i, x_j, ..., x_k) with cost values
	// E00...00, E00...01, E00...10,..., E11...10, E11...11.
	// Note the order of bits very carefully. It is defined consistent with the order
	// used in the QPBO software. E00...01 is the energy when only the LAST variable,
	// not the first, is 1.
	// vars is of size 'cliquesize'; E is of size 2^cliquesize
	// i < j < ... < k must hold
	void AddHigherTerm(int cliquesize, VID vars[], REAL E[])
	{
		int tn = 1 << cliquesize;
		for (int ix = 0; ix < tn; ix++)
		{
			REAL c = getCoef(ix, tn, E);
			if (c == 0)
				continue;
			VID vs[MAXD];
			int j = 0, b = 1 << (cliquesize - 1);
			for (int i = cliquesize - 1; i >= 0; i--)
			{
				if (ix & b)
					vs[j++] = vars[cliquesize - 1 - i];
				b >>= 1;
			}
			add(c, j, vs);
		}
	}

	// Delete terms with zero coefficient
	void shrink() {terms.shrink();}

	// Returns the maximum degree that MIGHT have any term.
	int maxD() const {return terms.maxD();}

	// Returns the maximum variable ID actually used.
	VID maxID() const {return terms.maxID();}

	// Reduces this PBF into a qpbf, adding variables (Tin's method). newvar is the variable ID for use as new variable.
	// Returns the variable ID to use for the next new variable.
	VID toQuadratic_Tin(PBF& qpbf, VID newvar, int W)
	{
		bool visited[newvar];
		for(int i = 0; i < newvar; i++)
			visited[i] = false;

		newvar = std::max(newvar, maxID()+1);
		for (iterator it = begin(); it != end(); ++it)
		{
			// Converts higher-order terms to quadratic (Tin's method), if necessary.
			int d = it.degree();
			variables va = it;
			VVecIt vi = va.begin();
			REAL c = it.coef();
			if (d == 0)
				qpbf.addConst(c);
			else if (d == 1)
				qpbf.add1(c, *vi);
			else if (d == 2)
				qpbf.add2(c, *vi, vi[1]);
			else if (d == 3)
			{	// Deal with d == 3 terms that are not found in the above cliques
				int idx = (*vi == vi[1] - W + 1) ? *vi - 1 : *vi;
				if(!visited[idx]){
					visited[idx] = true;

					VVec vars;
					vars.push_back(idx);
					vars.push_back(idx + 1);
					vars.push_back(idx + W);
					vars.push_back(idx + W + 1);
					VVecIt new_vi = vars.begin();

					smallPBF<REAL> clique;
					getClique(clique, variables(4, new_vi));
					addQuadraticTin(qpbf, clique, new_vi, newvar);
				}
			}
			else
			{
				visited[*vi] = true;

				smallPBF<REAL> clique;
				getClique(clique, it);
				addQuadraticTin(qpbf, clique, vi, newvar);
			}
		}
		qpbf.addConst(constant);
		return newvar;
	}

	void addQuadraticTin(PBF<REAL> &qpbf, smallPBF<REAL> &clique, VVecIt &vi, VID &newvar){
		bool flipped[4];
		for(int i = 0; i < 4; i++)	flipped[i] = false;

		std::vector<std::pair<BITMASK, REAL> > cubics;
		REAL quartic = 0;
		clique.get_hot_terms(cubics, quartic);

		smallPBF<REAL> extra;
		int lemma = apply_bitflips(clique, cubics, quartic, vi, extra, flipped);

		apply_lemma(lemma, cubics, quartic, vi, newvar, qpbf, flipped);
		qpbf.addQuad(extra, flipped, vi);
		/*
		VID v_true[4]; int perm[4];
		std::vector<std::pair<BITMASK, REAL> > cub;
		REAL quar = 0;
		clique.get_hot_terms(cub, quar);
		get_true_vi(cub, vi, v_true, perm);
		bool flag = check(quad, flipped, cub, quar, v_true, newvar);

		if(!flag){
			std::printf("%d, %d, %d, %d, quartic = %d\n", cubics[0].second, cubics[1].second, cubics[2].second, cubics[3].second, quartic);
//			std::printf("%d, %d, %d, %d\n", cubics[0].first, cubics[1].first, cubics[2].first, cubics[3].first);
			extra.print();
			system("PAUSE");
		}
		*/
		newvar++;
	}

	void apply_lemma(int lemma, std::vector<std::pair<BITMASK, REAL> > &cubics, REAL &quartic, VVecIt &vi, VID newvar, PBF<REAL> &quad, bool *flipped){
		VID vi_true[4];	int perm[4];
		get_true_vi(cubics, vi, vi_true, perm);
//		assert(check_lemma_condition(cubics, quartic, lemma));

		REAL cubics_sum = cubics[0].second + cubics[1].second + cubics[2].second + cubics[3].second;
		if(lemma == 1){
			quad.add1(3*quartic + cubics_sum, newvar);
			for (int i = 0; i < 3; i++)
				for (int j = i + 1; j < 4; j++)
					if((i == 0 && j == 3) || (i == 1 && j == 2))
						quad.add2_flipped(quartic + cubics_sum - cubics[i].second - cubics[j].second, vi_true[i], vi_true[j], flipped[perm[i]], flipped[perm[j]]);
					else
						quad.add2_flipped(quartic + cubics[i].second + cubics[j].second, vi_true[i], vi_true[j], flipped[perm[i]], flipped[perm[j]]);

			for(int i = 0; i < 4; i++)
				quad.add2_flipped( - 2*quartic - cubics_sum + cubics[3 - i].second, vi_true[i], newvar, flipped[perm[i]], false);
		}
		else{//Lemma 2
			quad.add1( - 3*quartic - 2*cubics_sum, newvar);
			for(int i = 0; i < 4; i++)
				quad.add2_flipped( quartic + cubics_sum - cubics[3 - i].second, vi_true[i], newvar, flipped[perm[i]], false);
		}
	}

	short int apply_bitflips(smallPBF<REAL> &clique,std::vector<std::pair<BITMASK, REAL> > &cubics, REAL &quartic, VVecIt &vi, smallPBF<REAL> &quad, bool *flipped){

		if(quartic < 0){
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			std::sort(cubics.begin(), cubics.end(), [](const std::pair<BITMASK, REAL> &a, const std::pair<BITMASK, REAL> &b) {return a.second < b.second;} );
		}

		unsigned short int case_ID = clique.get_Tin_case(cubics, quartic);	assert(case_ID != 0);
		switch(case_ID){
		case 3:
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 40:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 5:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 6:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 7:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 8:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 9:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 10:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 11:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 12:
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 14:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 160:
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 170:
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 180:
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 161:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 171:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 181:
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 19:
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		case 20:
			bitflip(cubics, quartic, quad, vi, flipped, 1);
			bitflip(cubics, quartic, quad, vi, flipped, 2);
			bitflip(cubics, quartic, quad, vi, flipped, 3);
			bitflip(cubics, quartic, quad, vi, flipped, 4);
			break;
		}

		if(case_ID == 3 || case_ID == 12)
			return 2;
		else
			return 1;
	}

	void add1_flipped(REAL c, VID v1, bool flip){
		if(flip){
			add1(-c, v1);
			addConst(c);
		}
		else
			add1(c, v1);
	}
	void add2_flipped(REAL c, VID v1, VID v2, bool flip1, bool flip2){
		if(v1 > v2){
			VID temp = v1;
			v1 = v2;
			v2 = temp;
			bool b = flip1;
			flip1 = flip2;
			flip2 = b;
		}

		if(flip1){
			if(flip2){
				add2(c, v1, v2);
				add1(-c, v1);
				add1(-c, v2);
				addConst(c);
			}
			else{
				add2(-c, v1, v2);
				add1(c, v2);
			}
		}
		else{
			if(flip2){
				add2(-c, v1, v2);
				add1(c, v1);
			}
			else
				add2(c, v1, v2);
		}
	}

	void addQuad(smallPBF<REAL> &extra, bool *flipped, VVecIt &vi){
		int v1, v2;
		for(int i = 0; i < 16; i++){
			v1 = i>>2;
			v2 = i%4;
			if(v1 > v2 || !extra.term[i]) continue;

			if(v1 == v2)
				add1_flipped(extra.term[i], vi[v1], flipped[v1]);
			else
				add2_flipped(extra.term[i], vi[v1], vi[v2], flipped[v1], flipped[v2]);
		}
		addConst(extra.cnst());
	}

	/*
	bool check_lemma_condition(std::vector<std::pair<BITMASK, REAL> > &cubics, REAL quartic, int lemma){
		if(lemma == 1){
			for(int i = 0; i < 4; i++){
				if(cubics[i].second < -quartic)
					return false;
				for(int j = i+1; j < 4; j++)
					if(cubics[i].second + cubics[j].second < -quartic)
						return false;
			}
			return true;
		}
		else if(lemma == 2){
			for(int i = 0; i < 4; i++)
				if(cubics[i].second > 0)
					return false;
			if(quartic > 0)
				return false;
			return true;
		}
	}

	bool check(PBF<REAL> &quad, bool *flipped, std::vector<std::pair<BITMASK, REAL> > &cubics, REAL quartic, VID *vi, VID newvar){
		bool flag = true;
		for(int i = 0; i < 16; i++){
			int bit[4];
			bit[0] = i%2;
			bit[1] = (i>>1)%2;
			bit[2] = (i>>2)%2;
			bit[3] = (i>>3)%2;

			REAL val0 = evaluate(quad, bit, vi, newvar, 0);
			REAL val1 = evaluate(quad, bit, vi, newvar, 1);
			REAL res = std::min(val0, val1);

			REAL org = 0, cubics_sum = cubics[0].second + cubics[1].second + cubics[2].second + cubics[3].second;
//			for(int k = 0; k < 4; k++)
//				if(flipped[vi[k]])
//					bit[k] = (bit[k] + 1) % 2;

			int k = bit[0] + 2*bit[1] + 4*bit[2] + 8*bit[3];
			if(k == 15)
				org = quartic + cubics_sum;
			else if(k == 14 || k == 13 || k == 11 || k == 7){
				for(int j = 0; j < 4; j++)
					if(bit[j] == 0)
						org = cubics[3-j].second;
			}
//			std::printf("val0 = %4d\tval1 = %4d\tmin = %4d\torg = %4d\n", val0, val1, res, org);
			if(org != res)
				flag = false;
		}
		return flag;
	}

	REAL evaluate(PBF<REAL> &quad, int *bit, VID *v, VID newvar, int aux){
		REAL value = 0;
//		std::printf("%d, %d, %d, %d\n", bit[0], bit[1], bit[2], bit[3]);
		for(iterator it = quad.begin(); it != quad.end(); ++it){
			int d = it.degree();
			variables va = it;
			VVecIt vi = va.begin();

			bool flag = true;
			for(int i = 0; i < d; i++){
//				std::printf("%d, ", vi[i]);
				for(int j = 0; j < 4; j++)
					if(vi[i] == v[j] && bit[j] == 0){
						flag = false;
						break;
					}
				if(aux == 0 && vi[i] == newvar)
					flag = false;
			}
			if(flag){
				value += it.coef();
//				std::printf("   Added %d", c);
			}
//			std::printf("\n");
		}
		value += quad.cnst();
//		std::printf("value = %d\n", value);
		return value;
	}

	*/
	void bitflip(std::vector< std::pair <BITMASK, REAL> > &cubics, REAL &quartic, smallPBF<REAL>& extra, VVecIt &vi_org, bool *flipped, short int bit){
		bit--;

		VID vi[4];
		int perm[4];
		get_true_vi(cubics, vi_org, vi, perm);
		/*
		std::printf("bitflip: %d\n", bit);
		std::printf("cubics: %d, %d, %d, %d\n", cubics[0].second, cubics[1].second, cubics[2].second, cubics[3].second);
		std::printf("perm: %d, %d, %d, %d\n", perm[0], perm[1], perm[2], perm[3]);
		*/
		flipped[perm[bit]] = !flipped[perm[bit]];

		int v1 = perm[bit];
		for(int i = v1<<2; i < (v1+1)<<2; i++){
			if(extra.term[i]){
				REAL c = extra.term[i];
				int v2 = i%4;

//				std::printf("v1 = %d, v2 = %d, c = %d\n", v1, v2, c);
				if(v1 == v2){
					extra.addConst(c);
					extra.term[i] = -c;
				}
				else{
					extra.term[(v2<<2) + v2] += c;
					assert(extra.term[(v2<<2) + v1] == c);
					extra.term[i] = -c;
					extra.term[(v2<<2) + v1] = -c;
				}
			}
		}

		if(bit == 0){
			extra.add2(cubics[0].second, perm[1], perm[2]);
			extra.add2(cubics[1].second, perm[1], perm[3]);
			extra.add2(cubics[2].second, perm[2], perm[3]);
		}
		else if(bit == 1){
			extra.add2(cubics[0].second, perm[0], perm[2]);
			extra.add2(cubics[1].second, perm[0], perm[3]);
			extra.add2(cubics[3].second, perm[2], perm[3]);
		}
		else if(bit == 2){
			extra.add2(cubics[0].second, perm[0], perm[1]);
			extra.add2(cubics[2].second, perm[0], perm[3]);
			extra.add2(cubics[3].second, perm[1], perm[3]);
		}
		else{// (bit == 3)
			extra.add2(cubics[1].second, perm[0], perm[1]);
			extra.add2(cubics[2].second, perm[0], perm[2]);
			extra.add2(cubics[3].second, perm[1], perm[2]);
		}


		for(int i = 0; i < 4; i++)
			if(bit == (3-i))
				cubics[i].second += quartic;
			else
				cubics[i].second *= -1;
		quartic *= -1;
		/*
		extra.print2();
		std::printf("cubics: %d, %d, %d, %d, quartic: %d\n", cubics[0].second, cubics[1].second, cubics[2].second, cubics[3].second, quartic);
		std::printf("    \nvi_org: %d, %d, %d, %d\n", vi_org[0], vi_org[1], vi_org[2], vi_org[3]);
		std::printf("bit: %d, vi: %d, %d, %d, %d\n", bit, vi[0], vi[1], vi[2], vi[3]);
		std::printf("%d, %d, %d, %d\n", cubics[0].first, cubics[1].first, cubics[2].first, cubics[3].first);
		system("PAUSE");
		*/
	}

	void get_true_vi(std::vector<std::pair<BITMASK, REAL> > &cubics, VVecIt &vi_org, VID *vi, int *perm){
		for(int i = 0; i < 4; i++)
			switch(cubics[i].first){
			case 7:
				vi[3-i] = vi_org[0];
				perm[3-i] = 0;
				break;
			case 11:
				vi[3-i] = vi_org[1];
				perm[3-i] = 1;
				break;
			case 13:
				vi[3-i] = vi_org[2];
				perm[3-i] = 2;
				break;
			case 14:
				vi[3-i] = vi_org[3];
				perm[3-i] = 3;
				break;
			default:
				std::printf("cubics %d bitmask not 7, 11, 13 or 14\n", i);
				system("PAUSE");
			}
//			cubics[3-i].second = cubics_org[3].second;

	}

	// Reduces this PBF into a qpbf, adding variables (HOCR). newvar is the variable ID for use as new variable.
	// Returns the variable ID to use for the next new variable.
	VID toQuadratic(PBF& qpbf, VID newvar) const
	{
		newvar = std::max(newvar, maxID()+1);
		for (iterator it = begin(); it != end(); ++it)
			qpbf.addQuadraticTerms(it, newvar);
		qpbf.addConst(constant);
		return newvar;
	}

	// Reduces higher-order terms (those with degree > 2) to quadratic by finding ELCs.
	// *Not all terms are always reduced.
	// When mindeg = 2 is specified, it additionally removes non-submodular quadratic terms as much as possible.
	void reduceHigher(int mindeg = 3)
	{
		for (iterator it = begin(); it != end(); ++it)
		{
			int d = it.degree();
			if (d < mindeg)
				break;
			REAL c = it.coef();
			if (d == 2 && c <= 0)
				continue;
			smallPBF<REAL> pbf;
			if (getStar(pbf, it))
			{
				bool parity = (d & 1) ? (c < 0) : (c > 0);
				REAL q = abs(c);
				int elc_idx = pbf.getELC(it, parity, q);
				if (elc_idx >= 0)
					addCancelTerm(it, elc_idx, q);
			}
		}
		shrink();
	}

	// Reduces higher-order terms (those with degree > 2) to quadratic by guessing ELCs.
	// All terms are always reduced, but *not all reductions are always correct*.
	void reduceHigherApprox()
	{
		for (iterator it = begin(); it != end(); ++it)
		{
			int d = it.degree();
			if (d < 3)
				break;
			REAL c = it.coef();
			bool parity = (d & 1) ? (c < 0) : (c > 0);
			smallPBF<REAL> pbf;
			getCore(pbf, it);
			int elc_idx = pbf.getMax(parity);
			addCancelTerm(it, elc_idx, abs(c));
		}
		shrink();
	}

	/*
	// enumerate terms

		for (iterator it = begin(); it != end(); ++it)
		{
			int degree = it.degree(); // the degree of the term
			VVecIt vars = it.vars(); // the variables
			REAL c = it.coef(); // the coefficient of the term
		}
	*/
	typedef typename Terms<REAL>::iterator iterator;
	iterator begin() const {return terms.begin();}
	iterator end() const {return terms.end();}

	// Returns the constant term.
	REAL cnst() const {return constant;}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Converts to an optimization object with the same interface as the QPBO object by
//	V. Kolmogorov
// OPTIMIZER: An optimizer class. The required interfaces are
// AddNode, AddUnaryTerm, and AddPairwiseTerm. See QPBO.h
// Provide varcount as the minimum variable ID that needs to be added to the QPBO object.
// (i.e., if you are going to inquire the QPBO object if the variable has 0 or 1.)

	template<typename OPTIMIZER>
	void convert(OPTIMIZER& optimization, int varcount, FILE* Log)
	{
		int mini = 100000, maxi = -100000, quadratics = 0;

		varcount = std::max(varcount, maxID() + 1);
		optimization.AddNode(varcount);
		for (iterator it = begin(); it != end(); ++it)
		{
			int d = it.degree();
			if (d > 2)
				exit(-1);
			VVecIt vs = it.vars();
			REAL c = it.coef();

			if (c < mini) mini = c;
			if (c > maxi) maxi = c;

			if (d == 1)
				optimization.AddUnaryTerm(*vs, 0, c);
			else{
				optimization.AddPairwiseTerm(*vs, vs[1], 0, 0, 0, c);
				quadratics++;
			}
		}
		optimization.AddUnaryTerm(0, cnst(), cnst());

		printf("\t[%7d, %7d]\t%d\t\t", mini, maxi, quadratics);
		fprintf(Log,"\t[%7d, %7d]\t%d\t\t", mini, maxi, quadratics);
	}

private:
	// Adds a monomial.
	void addConst(REAL c) {constant += c;}
	void add(REAL c, int degree, VID* vars) {if (degree == 0) constant += c; else terms.add(c, degree, vars);}
	void add1(REAL c, VID v) {terms.add1(c, v);}
	void add2(REAL c, VID v1, VID v2) {terms.add2(c, v1, v2);}


	// Converts higher-order terms to quadratic (HOCR), if necessary.
	void addQuadraticTerms(const iterator& it, VID& newvar)
	{
		int d = it.degree();
		variables va = it;
		VVecIt vi = va.begin();
		REAL c = it.coef();
		if (d == 0)
			addConst(c);
		else if (d == 1)
			add1(c, *vi);
		else if (d == 2)
			add2(c, *vi, vi[1]);
		else if (c < 0)
		{
			for (; vi != va.end(); ++vi)
				add2(c, *vi, newvar);
			add1((1 - d) * c, newvar);
			newvar++;
		}
		else
		{
			int numNewVars = ((d % 2) ? (d - 1) : (d - 2)) / 2;
			for (int i = 1; i <= numNewVars; i++)
			{
				bool b = (d % 2) && i == numNewVars;
				REAL coef = c * (b ? -1 : -2);
				for (int j = 0; j < d; j++) // S_1
					add2(coef, vi[j], newvar);
				add1(c * ((b ? 2 : 4) * i - 1), newvar);
				newvar++;
			}
			for (int i = 0; i < d - 1; i++) // S_2
				for (int j = i + 1; j < d; j++)
					add2(c, vi[i], vi[j]);
		}
	}
	// cancel the highest order term.
	void addCancelTerm(const variables& vars, int elc_idx, REAL C)
	{
		int cliquesize = vars.size();
		int tn = 1 << cliquesize;
		for (int ix = 0; ix < tn; ix++)
		{
			REAL c = C * matrixElement(ix, elc_idx);
			if (c == 0)
				continue;
			VID vs[MAXD];
			int j = 0, b = 1 << (cliquesize - 1);
			for (VVecIt vi = vars.begin(); vi != vars.end(); ++vi)
			{
				if (ix & b)
					vs[j++] = *vi;
				b >>= 1;
			}
			add(c, j, vs);
		}
	}
	// get all the terms that contains at least one of the variables given in vars
	bool getStar(smallPBF<REAL>& rv, const variables& core)
	{
		std::vector<TermPointer> ps;
		for (VVecIt vi = core.begin(); vi != core.end(); ++vi)
		{
			std::vector<TermPointer>& ts = terms.getTermsContaining(*vi);
			ps.insert(ps.end(), ts.begin(), ts.end());
		}
		std::sort(ps.begin(), ps.end());
		auto it2 = std::unique(ps.begin(), ps.end());
		std::vector<VID> vs;
		core.copy(std::back_inserter(vs));
		for (auto it = ps.begin(); it != it2; it++)
			if (terms.coef(*it))
			{
				variables va(terms, *it);
				va.copy(std::back_inserter(vs));
			}
		std::sort(vs.begin(), vs.end());
		auto it3 = std::unique(vs.begin(), vs.end());
		int numvars = (int)(it3 - vs.begin());
		if (numvars > MAXVARELC + core.size())
			return false; // too big
		rv.setVars(numvars, vs.begin());
		for (auto it = ps.begin(); it != it2; it++)
			if (terms.coef(*it))
				rv.add(variables(terms, *it), terms.coef(*it));
		rv.addConst(constant);
		return true;
	}

	// get all the terms that consist of only the variables given in vars
	void getCore(smallPBF<REAL>& rv, const variables& vars)
	{
		rv.setVars(vars);
		std::vector<TermPointer> ps;
		for (VVecIt vi = vars.begin(); vi != vars.end(); ++vi)
		{
			std::vector<TermPointer>& ts = terms.getTermsContaining(*vi);
			for (auto it = ts.begin(); it != ts.end(); it++)
				if (vars.includes(variables(terms, *it)))
					ps.push_back(*it);
		}
		std::sort(ps.begin(), ps.end());
		auto it2 = std::unique(ps.begin(), ps.end());
		for (auto it = ps.begin(); it != it2; it++)
			rv.add(variables(terms, *it), terms.coef(*it));
		rv.addConst(constant);
	}

	void getClique(smallPBF<REAL>& rv, const variables& vars)
	{
		rv.setVars(vars);
		std::vector<TermPointer> ps;
		for (VVecIt vi = vars.begin(); vi != vars.end(); ++vi)
		{
			std::vector<TermPointer>& ts = terms.getTermsContaining(*vi);
			for (auto it = ts.begin(); it != ts.end(); it++){
				if(it->degree() < 3)	continue;
				if (vars.includes(variables(terms, *it)))
					ps.push_back(*it);
			}
		}
		std::sort(ps.begin(), ps.end());
		auto it2 = std::unique(ps.begin(), ps.end());
		for (auto it = ps.begin(); it != it2; it++)
			rv.add(variables(terms, *it), terms.coef(*it));
	}

	REAL constant;
	Terms<REAL> terms;
};




} // namespace ELCReduce


#endif // __ELC_H__

