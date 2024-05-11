// Copyright (C) 2015 Chris Sweeney. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Chris Sweeney nor the names of its contributors may
//       be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include <trajectory_utils/polynomial.h>

#include <Eigen/Core>
#include <cmath>
#include <iostream>

using namespace std;


namespace rpoly_plus_plus {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

// Remove leading terms with zero coefficients.
VectorXd RemoveLeadingZeros(const VectorXd& polynomial_in) {
  int i = 0;
  while (i < (polynomial_in.size() - 1) && polynomial_in(i) == 0) {
    ++i;
  }
  return polynomial_in.tail(polynomial_in.size() - i);
}

VectorXd DifferentiatePolynomial(const VectorXd& polynomial) {
  const int degree = polynomial.rows() - 1;

  // Degree zero polynomials are constants, and their derivative does
  // not result in a smaller degree polynomial, just a degree zero
  // polynomial with value zero.
  if (degree == 0) {
    return VectorXd::Zero(1);
  }

  VectorXd derivative(degree);
  for (int i = 0; i < degree; ++i) {
    derivative(i) = (degree - i) * polynomial(i);
  }

  return derivative;
}

VectorXd MultiplyPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  VectorXd multiplied_poly = VectorXd::Zero(poly1.size() + poly2.size() - 1);;
  for (int i = 0; i < poly1.size(); i++) {
    for (int j = 0; j < poly2.size(); j++) {
      multiplied_poly.reverse()(i + j) +=
          poly1.reverse()(i) * poly2.reverse()(j);
    }
  }
  return multiplied_poly;
}

VectorXd AddPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  if (poly1.size() > poly2.size()) {
    VectorXd sum = poly1;
    sum.tail(poly2.size()) += poly2;
    return sum;
  } else {
    VectorXd sum = poly2;
    sum.tail(poly1.size()) += poly1;
    return sum;
  }
}

double FindRootIterativeNewton(const Eigen::VectorXd& polynomial,
                               const double x0,
                               const double epsilon,
                               const int max_iterations) {
  double root = x0;
  const Eigen::VectorXd derivative = DifferentiatePolynomial(polynomial);
  double prev;
  for (int i = 0; i < max_iterations; i++) {
    prev = root;
    root -= EvaluatePolynomial(polynomial, root) /
            EvaluatePolynomial(derivative, root);
    if (std::abs(prev - root) < epsilon) {
      break;
    }
  }
  return root;
}

// Find the remainder of a poly1/poly2.
Eigen::VectorXd RemPolynomials(const Eigen::VectorXd& poly1,
                               const Eigen::VectorXd& poly2){
	//If polynomial 1 is less than polynomial 2 just return
	Eigen::VectorXd divider = poly1;
	Eigen::VectorXd remainder = divider;
	//std::cout << "poly2 " << poly2.size() <<std::endl;
	//std::cout << "divider " << divider.size() <<std::endl;
	//Keep dividing till it can no longer be divided by this 
	while (divider.size() >= poly2.size()){
		double factor = divider[0]/poly2[0];
		remainder = VectorXd::Zero(divider.size()-1);
		for(int i= 1; i < poly2.size();i++){
			remainder[i-1] = divider[i] - factor*poly2[i];
		}	
		//Fill in the remainder 
		for(int i = poly2.size(); i < divider.size(); i++){
			remainder[i-1] = divider[i];
		}			
		remainder = RemoveLeadingZeros(remainder);
		divider = remainder; 
	}	
	return remainder;
}

int FindSturmRoot(const Eigen::VectorXd& poly1,
				   const double start,
                   const double end){
	int initSeq = 0;
	int endSeq = 0;
	double init2Val;
	double end2Val;
	double init1Val = EvaluatePolynomial(poly1, start);
	double end1Val = EvaluatePolynomial(poly1, end);
	Eigen::VectorXd prevPoly = poly1;
	//std::cout << poly1 << std::endl;
	Eigen::VectorXd currPoly = DifferentiatePolynomial(poly1);
	while(currPoly.size()>1){
	//for(int i =0;i<2;i++){
		//std::cout << "init1VAL " << init1Val << std::endl;
		//std::cout << "end1Val " << end1Val << std::endl;
		init2Val = EvaluatePolynomial(currPoly, start);
		end2Val = EvaluatePolynomial(currPoly, end);
		//Detect sign flip
		if(init2Val*init1Val < 0){
			initSeq+=1;
		}
		if(end2Val*end1Val < 0){
			endSeq+=1;
		}
		//Copy paste new value
		if(init2Val !=0){
			init1Val = init2Val;
		}
		if(end2Val !=0){
			end1Val = end2Val;
		}
		//std::cout << "poly1 " << currPoly << std::endl;
		//std::cout << "divider " <<prevPoly<< std::endl;
		Eigen::VectorXd currPoly2 = -1*RemPolynomials(prevPoly, currPoly  );
		prevPoly = currPoly;
		currPoly = currPoly2;
	//}
	}
	if(init1Val*currPoly[0] < 0){
		initSeq+=1;
	}
	if(end1Val*currPoly[0] < 0){
		endSeq+=1;
	}
	return initSeq - endSeq;
}

}  // namespace rpoly_plus_plus
