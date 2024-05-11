/*
    MIT License

    Copyright (c) 2020 Zhepei Wang (wangzhepei@live.com)
                  2021 Hongkai Ye (kyle_yeh@163.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
#pragma once

#include <Eigen/Eigen>
#include <cmath>
#include <vector>

#include "root_finder.hpp"

// Polynomial order and trajectory dimension are fixed here
constexpr int TrajOrder = 5;
constexpr int TrajDim = 3;

// Type for piece boundary condition and coefficient matrix
typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> BoundaryCond;
typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> CoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder> VelCoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder - 1> AccCoefficientMat;
typedef Eigen::Matrix<double, TrajDim, TrajOrder - 2> JerkCoefficientMat;

typedef Eigen::Matrix<double, 6, 1> StatePV;
typedef Eigen::Matrix<double, 9, 1> StatePVA;
typedef Eigen::Matrix<double, TrajDim, 1> ControlJrk;
typedef Eigen::Matrix<double, TrajDim, 1> ControlAcc;

typedef Eigen::Matrix<double, 1, TrajOrder + 1> YawBoundaryCond;
typedef Eigen::Matrix<double, 1, TrajOrder + 1> YawCoefficientMat;
typedef Eigen::Matrix<double, 1, TrajOrder> YawDotCoefficientMat;
typedef Eigen::Matrix<double, 1, TrajOrder - 1> YawDdotCoefficientMat;

// A single piece of a trajectory, which is indeed a polynomial
class Piece
{
private:
  // Piece(t) = c5*t^5 + c4*t^4 + ... + c1*t + c0
  // The natural coefficient matrix = [c5,c4,c3,c2,c1,c0]
  double duration;
  CoefficientMat CoeffMat;

public:
  Piece() = default;

  // Constructor from duration and coefficient
  Piece(double dur, CoefficientMat coeffs) : duration(dur), CoeffMat(coeffs)
  {
  }

  // Constructor from boundary condition and duration
  Piece(BoundaryCond boundCond, double dur) : duration(dur)
  {
    // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
    double t1 = dur;
    double t2 = t1 * t1;
    CoefficientMat nCoeffMat;

    // Inverse mapping is computed without explicit matrix inverse
    // It maps boundary condition to normalized coefficient matrix
    nCoeffMat.col(0) = 0.5 * (boundCond.col(5) - boundCond.col(2)) * t2 -
                       3.0 * (boundCond.col(1) + boundCond.col(4)) * t1 +
                       6.0 * (boundCond.col(3) - boundCond.col(0));
    nCoeffMat.col(1) = (-boundCond.col(5) + 1.5 * boundCond.col(2)) * t2 +
                       (8.0 * boundCond.col(1) + 7.0 * boundCond.col(4)) * t1 +
                       15.0 * (-boundCond.col(3) + boundCond.col(0));
    nCoeffMat.col(2) = (0.5 * boundCond.col(5) - 1.5 * boundCond.col(2)) * t2 -
                       (6.0 * boundCond.col(1) + 4.0 * boundCond.col(4)) * t1 +
                       10.0 * (boundCond.col(3) - boundCond.col(0));
    nCoeffMat.col(3) = 0.5 * boundCond.col(2) * t2;
    nCoeffMat.col(4) = boundCond.col(1) * t1;
    nCoeffMat.col(5) = boundCond.col(0);

    double t = 1.0;
    for (int i = TrajOrder; i >= 0; i--)
    {
      CoeffMat.col(i) = nCoeffMat.col(i) / t;
      t *= dur;
    }
  }

  inline int getDim() const
  {
    return TrajDim;
  }

  inline int getOrder() const
  {
    return TrajOrder;
  }

  inline double getDuration() const
  {
    return duration;
  }

  inline void setDuration(double duration)
  {
    this->duration = duration;
  }

  // Get the position at time t in this piece
  inline Eigen::Vector3d getPos(double t) const
  {
    // Normalize the time
    Eigen::Vector3d pos(0.0, 0.0, 0.0);
    double tn(1.0);
    for (int i = TrajOrder; i >= 0; i--)
    {
      pos += tn * CoeffMat.col(i);
      tn *= t;
    }
    // The pos is not affected by normalization
    return pos;
  }

  // Get the velocity at time t in this piece
  inline Eigen::Vector3d getVel(double t) const
  {
    // Normalize the time
    Eigen::Vector3d vel(0.0, 0.0, 0.0);
    double tn(1.0);
    int n = 1;
    for (int i = TrajOrder - 1; i >= 0; i--)
    {
      vel += n * tn * CoeffMat.col(i);
      tn *= t;
      n++;
    }
    // Recover the actual vel
    return vel;
  }

  // Get the acceleration at time t in this piece
  inline Eigen::Vector3d getAcc(double t) const
  {
    // Normalize the time
    Eigen::Vector3d acc(0.0, 0.0, 0.0);
    double tn(1.0);
    int m = 1;
    int n = 2;
    for (int i = TrajOrder - 2; i >= 0; i--)
    {
      acc += m * n * tn * CoeffMat.col(i);
      tn *= t;
      m++;
      n++;
    }
    // Recover the actual acc
    return acc;
  }

  // Get the jerk at time t in this piece
  inline Eigen::Vector3d getJerk(double t) const
  {
    // Normalize the time
    Eigen::Vector3d jerk(0.0, 0.0, 0.0);
    double tn(1.0);
    int m = 1;
    int n = 2;
    int k = 3;
    for (int i = TrajOrder - 3; i >= 0; i--)
    {
      jerk += k * m * n * tn * CoeffMat.col(i);
      tn *= t;
      k++;
      m++;
      n++;
    }
    // Recover the actual acc
    return jerk;
  }

  // Get the snap at time t in this piece
  inline Eigen::Vector3d getSnap(double t) const
  {
    // Normalize the time
    Eigen::Vector3d snap(0.0, 0.0, 0.0);
    double tn(1.0);
    int m = 1;
    int n = 2;
    int k = 3;
    int w = 4;
    for (int i = TrajOrder - 4; i >= 0; i--)
    {
      snap += w * k * m * n * tn * CoeffMat.col(i);
      tn *= t;
      w++;
      k++;
      m++;
      n++;
    }
    // Recover the actual acc
    return snap;
  }

  // Get the boundary condition of this piece
  inline BoundaryCond getBoundCond() const
  {
    BoundaryCond boundCond;
    boundCond << getPos(0.0), getVel(0.0), getAcc(0.0),
        getPos(duration), getVel(duration), getAcc(duration);
    return boundCond;
  }

  // Get the coefficient matrix of the piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline CoefficientMat getCoeffMat(bool normalized = false) const
  {
    CoefficientMat posCoeffsMat;
    double t = 1;
    for (int i = TrajOrder; i >= 0; i--)
    {
      posCoeffsMat.col(i) = CoeffMat.col(i) * t;
      t *= normalized ? duration : 1.0;
    }
    return posCoeffsMat;
  }

  // Get the polynomial coefficients of velocity of this piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline VelCoefficientMat getVelCoeffMat(bool normalized = false) const
  {
    VelCoefficientMat velCoeffMat;
    int n = 1;
    double t = 1.0;
    t *= normalized ? duration : 1.0;
    for (int i = TrajOrder - 1; i >= 0; i--)
    {
      velCoeffMat.col(i) = n * CoeffMat.col(i) * t;
      n++;
      t *= normalized ? duration : 1.0;
    }
    return velCoeffMat;
  }

  // Get the polynomial coefficients of acceleration of this piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline AccCoefficientMat getAccCoeffMat(bool normalized = false) const
  {
    AccCoefficientMat accCoeffMat;
    int n = 2;
    int m = 1;
    double t = 1.0;
    t *= normalized ? duration * duration : 1.0;
    for (int i = TrajOrder - 2; i >= 0; i--)
    {
      accCoeffMat.col(i) = n * m * CoeffMat.col(i) * t;
      n++;
      m++;
      t *= normalized ? duration : 1.0;
    }
    return accCoeffMat;
  }

  // Get the polynomial coefficients of jerk of this piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline JerkCoefficientMat getJerkCoeffMat(bool normalized = false) const
  {
    JerkCoefficientMat jerkCoeffMat;
    int n = 3;
    int m = 2;
    int l = 1;
    double t = 1.0;
    t *= normalized ? duration * duration * duration : 1.0;
    for (int i = TrajOrder - 3; i >= 0; i--)
    {
      jerkCoeffMat.col(i) = n * m * l * CoeffMat.col(i) * t;
      n++;
      m++;
      l++;
      t *= normalized ? duration : 1.0;
    }
    return jerkCoeffMat;
  }

  // Get the max velocity rate of the piece
  inline double getMaxVelRate() const
  {
    // Compute normalized squared vel norm polynomial coefficient matrix
    Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
    Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                            RootFinder::polySqr(nVelCoeffMat.row(1)) +
                            RootFinder::polySqr(nVelCoeffMat.row(2));
    int N = coeff.size();
    int n = N - 1;
    for (int i = 0; i < N; i++)
    {
      coeff(i) *= n;
      n--;
    }
    if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
    {
      return 0.0;
    }
    else
    {
      // Search an open interval whose boundaries are not zeros
      double l = -0.0625;
      double r = 1.0625;
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
      {
        l = 0.5 * l;
      }
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
      {
        r = 0.5 * (r + 1.0);
      }
      // Find all stationaries
      std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                FLT_EPSILON / duration);

      // Check boundary points and stationaries within duration
      candidates.insert(0.0);
      candidates.insert(1.0);
      double maxVelRateSqr = -INFINITY;
      double tempNormSqr;
      for (std::set<double>::const_iterator it = candidates.begin();
           it != candidates.end();
           it++)
      {
        if (0.0 <= *it && 1.0 >= *it)
        {
          // Recover the actual time then get the vel squared norm
          tempNormSqr = getVel((*it) * duration).squaredNorm();
          maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
        }
      }
      return sqrt(maxVelRateSqr);
    }
  }

  // Get the max acceleration rate of the piece
  inline double getMaxAccRate() const
  {
    // Compute normalized squared acc norm polynomial coefficient matrix
    Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
    Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                            RootFinder::polySqr(nAccCoeffMat.row(1)) +
                            RootFinder::polySqr(nAccCoeffMat.row(2));
    int N = coeff.size();
    int n = N - 1;
    for (int i = 0; i < N; i++)
    {
      coeff(i) *= n;
      n--;
    }
    if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
    {
      return 0.0;
    }
    else
    {
      // Search an open interval whose boundaries are not zeros
      double l = -0.0625;
      double r = 1.0625;
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
      {
        l = 0.5 * l;
      }
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
      {
        r = 0.5 * (r + 1.0);
      }
      // Find all stationaries
      std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                FLT_EPSILON / duration);
      // Check boundary points and stationaries within duration
      candidates.insert(0.0);
      candidates.insert(1.0);
      double maxAccRateSqr = -INFINITY;
      double tempNormSqr;
      for (std::set<double>::const_iterator it = candidates.begin();
           it != candidates.end();
           it++)
      {
        if (0.0 <= *it && 1.0 >= *it)
        {
          // Recover the actual time then get the acc squared norm
          tempNormSqr = getAcc((*it) * duration).squaredNorm();
          maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
        }
      }
      return sqrt(maxAccRateSqr);
    }
  }

  inline double getMaxJerkRate() const
  {
    // Compute normalized squared jerk norm polynomial coefficient matrix
    Eigen::MatrixXd nJerkCoeffMat = getJerkCoeffMat(true);
    Eigen::VectorXd coeff = RootFinder::polySqr(nJerkCoeffMat.row(0)) +
                            RootFinder::polySqr(nJerkCoeffMat.row(1)) +
                            RootFinder::polySqr(nJerkCoeffMat.row(2));
    int N = coeff.size();
    int n = N - 1;
    for (int i = 0; i < N; i++)
    {
      coeff(i) *= n;
      n--;
    }
    if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
    {
      return 0.0;
    }
    else
    {
      // Search an open interval whose boundaries are not zeros
      double l = -0.0625;
      double r = 1.0625;
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
      {
        l = 0.5 * l;
      }
      while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
      {
        r = 0.5 * (r + 1.0);
      }
      // Find all stationaries
      std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                FLT_EPSILON / duration);
      // Check boundary points and stationaries within duration
      candidates.insert(0.0);
      candidates.insert(1.0);
      double maxJerkRateSqr = -INFINITY;
      double tempNormSqr;
      for (std::set<double>::const_iterator it = candidates.begin();
           it != candidates.end();
           it++)
      {
        if (0.0 <= *it && 1.0 >= *it)
        {
          // Recover the actual time then get the acc squared norm
          tempNormSqr = getJerk((*it) * duration).squaredNorm();
          maxJerkRateSqr = maxJerkRateSqr < tempNormSqr ? tempNormSqr : maxJerkRateSqr;
        }
      }
      return sqrt(maxJerkRateSqr);
    }
  }

  // Check whether velocity rate of the piece is always less than maxVelRate
  inline bool checkMaxVelRate(double maxVelRate) const
  {
    double sqrMaxVelRate = maxVelRate * maxVelRate;
    if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
        getVel(duration).squaredNorm() >= sqrMaxVelRate)
    {
      return false;
    }
    else
    {
      Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
      Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                              RootFinder::polySqr(nVelCoeffMat.row(1)) +
                              RootFinder::polySqr(nVelCoeffMat.row(2));
      // Convert the actual squared maxVelRate to a normalized one
      double t2 = duration * duration;
      coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
      // Directly check the root existence in the normalized interval
      return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
    }
  }

  // Check whether accleration rate of the piece is always less than maxAccRate
  inline bool checkMaxAccRate(double maxAccRate) const
  {
    double sqrMaxAccRate = maxAccRate * maxAccRate;
    if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
        getAcc(duration).squaredNorm() >= sqrMaxAccRate)
    {
      return false;
    }
    else
    {
      Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
      Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                              RootFinder::polySqr(nAccCoeffMat.row(1)) +
                              RootFinder::polySqr(nAccCoeffMat.row(2));
      // Convert the actual squared maxAccRate to a normalized one
      double t2 = duration * duration;
      double t4 = t2 * t2;
      coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
      // Directly check the root existence in the normalized interval
      return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
    }
  }

  // Check whether jerk rate of the piece is always less than maxJerkRate
  inline bool checkMaxJerkRate(double maxJerkRate) const
  {
    double sqrMaxJerkRate = maxJerkRate * maxJerkRate;
    if (getJerk(0.0).squaredNorm() >= sqrMaxJerkRate ||
        getJerk(duration).squaredNorm() >= sqrMaxJerkRate)
    {
      return false;
    }
    else
    {
      Eigen::MatrixXd nJerkCoeffMat = getJerkCoeffMat(true);
      Eigen::VectorXd coeff = RootFinder::polySqr(nJerkCoeffMat.row(0)) +
                              RootFinder::polySqr(nJerkCoeffMat.row(1)) +
                              RootFinder::polySqr(nJerkCoeffMat.row(2));
      // Convert the actual squared maxJerkRate to a normalized one
      double t2 = duration * duration;
      double t4 = t2 * t2;
      double t6 = t4 * t2;
      coeff.tail<1>()(0) -= sqrMaxJerkRate * t6;
      // Directly check the root existence in the normalized interval
      return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
    }
  }

  // Scale the Piece(t) to Piece(k*t)
  inline void scaleTime(double k)
  {
    duration /= k;
    return;
  }

  inline void sampleOneSeg(std::vector<StatePVA> *vis_x) const
  {
    double dt = 0.005;
    for (double t = 0.0; t < duration; t += dt)
    {
      Eigen::Vector3d pos, vel, acc;
      pos = getPos(t);
      vel = getVel(t);
      acc = getAcc(t);
      StatePVA x;
      x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
      vis_x->push_back(x);
    }
  }

  // for 5th degree polynomial
  inline void cutPiece(const Piece &orig_piece, double ts, CoefficientMat &new_coeff) const
  {
    CoefficientMat ori_coeff = orig_piece.getCoeffMat();
    double ts2 = ts * ts;
    double ts3 = ts2 * ts;
    double ts4 = ts3 * ts;
    double ts5 = ts4 * ts;
    for (int dim = 0; dim < 3; ++dim)
    {
      new_coeff(dim, 0) = ori_coeff(dim, 0);                              // c5*t^5
      new_coeff(dim, 1) = ori_coeff(dim, 1) + 5 * ori_coeff(dim, 0) * ts; // c4*4^4
      new_coeff(dim, 2) = ori_coeff(dim, 2) + 4 * ori_coeff(dim, 1) * ts + 10 * ori_coeff(dim, 0) * ts2;
      new_coeff(dim, 3) = ori_coeff(dim, 3) + 3 * ori_coeff(dim, 2) * ts + 6 * ori_coeff(dim, 1) * ts2 + 10 * ori_coeff(dim, 0) * ts3;
      new_coeff(dim, 4) = ori_coeff(dim, 4) + 2 * ori_coeff(dim, 3) * ts + 3 * ori_coeff(dim, 2) * ts2 + 4 * ori_coeff(dim, 1) * ts3 + 5 * ori_coeff(dim, 0) * ts4;
      new_coeff(dim, 5) = ori_coeff(dim, 5) + ori_coeff(dim, 4) * ts + ori_coeff(dim, 3) * ts2 + ori_coeff(dim, 2) * ts3 + ori_coeff(dim, 1) * ts4 + ori_coeff(dim, 0) * ts5;
    }
  }

  // for 5th degree polynomial
  inline double calCost(const double &rho) const
  {
    AccCoefficientMat acc_coe = getAccCoeffMat();
    double cost(0.0);
    int n_dim = acc_coe.rows();
    double tau2 = duration * duration;
    double tau3 = tau2 * duration;
    double tau4 = tau3 * duration;
    double tau5 = tau4 * duration;
    double tau6 = tau5 * duration;
    for (int dim = 0; dim < n_dim; ++dim)
    {
      double cost_curr_dim = acc_coe(dim, 0) * acc_coe(dim, 0) * tau6 / 7.0 + acc_coe(dim, 0) * acc_coe(dim, 1) * tau5 / 3.0 + acc_coe(dim, 0) * acc_coe(dim, 2) * tau4 * 2.0 / 5.0 + acc_coe(dim, 0) * acc_coe(dim, 3) * tau3 / 2.0 + acc_coe(dim, 1) * acc_coe(dim, 1) * tau4 / 5.0 + acc_coe(dim, 1) * acc_coe(dim, 2) * tau3 / 2.0 + acc_coe(dim, 1) * acc_coe(dim, 3) * tau2 * 2.0 / 3.0 + acc_coe(dim, 2) * acc_coe(dim, 2) * tau2 / 3.0 + acc_coe(dim, 2) * acc_coe(dim, 3) * duration + acc_coe(dim, 3) * acc_coe(dim, 3);
      cost += cost_curr_dim;
    }
    cost *= rho;
    cost += 1.0;
    cost *= duration;
    return cost;
  }

  inline double getTotalDistance(const double &dt) const
  {
    double dist(0.0);
    Eigen::Vector3d p0 = getPos(0.0), p;
    for (double t = dt; t <= duration; t += dt)
    {
      p = getPos(t);
      dist += (p - p0).norm();
      p0 = p;
    }
    return dist;
  }

  inline double getDistance(const double &dt, const double &s, const double &e) const
  {
    double dist(0.0);
    Eigen::Vector3d p0 = getPos(s), p;
    for (double t = s + dt; t <= e; t += dt)
    {
      p = getPos(t);
      dist += (p - p0).norm();
      p0 = p;
    }
    return dist;
  }

  inline double project_pt(const Eigen::Vector3d &pt,
                           double &tt, Eigen::Vector3d &pro_pt)
  {
    // 2*(p-p0)^T * \dot{p} = 0
    auto l_coeff = getCoeffMat();
    l_coeff.col(5) = l_coeff.col(5) - pt;
    auto r_coeff = getVelCoeffMat();
    Eigen::VectorXd eq = Eigen::VectorXd::Zero(2 * 5);
    for (int j = 0; j < l_coeff.rows(); ++j)
    {
      eq = eq + RootFinder::polyConv(l_coeff.row(j), r_coeff.row(j));
    }
    double l = -0.0625;
    double r = duration + 0.0625;
    while (fabs(RootFinder::polyVal(eq, l)) < DBL_EPSILON)
    {
      l = 0.5 * l;
    }
    while (fabs(RootFinder::polyVal(eq, r)) < DBL_EPSILON)
    {
      r = 0.5 * (duration + r);
    }
    std::set<double> roots =
        RootFinder::solvePolynomial(eq, l, r, 1e-6);
    // std::cout << "# roots: " << roots.size() << std::endl;
    double min_dist = -1;
    for (const auto &root : roots)
    {
      // std::cout << "root: " << root << std::endl;
      if (root < 0 || root > duration)
      {
        continue;
      }
      if (getVel(root).norm() < 1e-6)
      { // velocity == 0, ignore it
        continue;
      }
      Eigen::Vector3d p = getPos(root);
      double distance = (p - pt).norm();
      if (distance < min_dist || min_dist < 0)
      {
        min_dist = distance;
        tt = root;
        pro_pt = p;
      }
    }
    return min_dist;
  }

  inline std::vector<Eigen::Vector3d> get_all_project_pt(const Eigen::Vector3d &pt, const double &t_start, std::vector<double> &best_t_set)
  {
    // 2*(p-p0)^T * \dot{p} = 0
    std::vector<Eigen::Vector3d> pro_pt_set;
    best_t_set.clear();
    auto l_coeff = getCoeffMat();
    l_coeff.col(5) = l_coeff.col(5) - pt;
    auto r_coeff = getVelCoeffMat();
    Eigen::VectorXd eq = Eigen::VectorXd::Zero(2 * 5);
    for (int j = 0; j < l_coeff.rows(); ++j)
    {
      eq = eq + RootFinder::polyConv(l_coeff.row(j), r_coeff.row(j));
    }
    double l = -0.0625;
    double r = duration + 0.0625;
    while (fabs(RootFinder::polyVal(eq, l)) < DBL_EPSILON)
    {
      l = 0.5 * l;
    }
    while (fabs(RootFinder::polyVal(eq, r)) < DBL_EPSILON)
    {
      r = 0.5 * (duration + r);
    }
    std::set<double> roots =
        RootFinder::solvePolynomial(eq, l, r, 1e-6);
    for (const auto &root : roots)
    {
      // std::cout << "root: " << root << std::endl;
      if (root < t_start || root < 0 || root > duration)
      {
        continue;
      }
      if (getVel(root).norm() < 1e-6)
      { // velocity == 0, ignore it
        continue;
      }
      // std::cout << "find min!" << std::endl;
      Eigen::Vector3d p = getPos(root);
      Eigen::Vector3d pv = getVel(root);
      if ((p - pt).dot(pv) > 1e-4)
      {
        std::cout << "wrong root occur !" << std::endl;
        continue;
      }
      pro_pt_set.push_back(p);
      best_t_set.push_back(root);
      // std::cout << "p: " << p.transpose() << std::endl;
    }
    return pro_pt_set;
  }
};

class Piece1D
{
private:
  // Piece(t) = c5*t^5 + c4*t^4 + ... + c1*t + c0
  // The natural coefficient matrix = [c5,c4,c3,c2,c1,c0]
  double duration;
  YawCoefficientMat CoeffMat;

public:
  Piece1D() = default;

  // Constructor from duration and coefficient
  Piece1D(double dur, YawCoefficientMat coeffs) : duration(dur), CoeffMat(coeffs)
  {
  }

  // Constructor from boundary condition and duration
  Piece1D(YawBoundaryCond boundCond, double dur) : duration(dur)
  {
    // The BoundaryCond matrix boundCond = [p(0),v(0),a(0),p(T),v(T),a(T)]
    double t1 = dur;
    double t2 = t1 * t1;
    YawCoefficientMat nCoeffMat;

    // Inverse mapping is computed without explicit matrix inverse
    // It maps boundary condition to normalized coefficient matrix
    nCoeffMat(0) = 0.5 * (boundCond(5) - boundCond(2)) * t2 -
                   3.0 * (boundCond(1) + boundCond(4)) * t1 +
                   6.0 * (boundCond(3) - boundCond(0));
    nCoeffMat(1) = (-boundCond(5) + 1.5 * boundCond(2)) * t2 +
                   (8.0 * boundCond(1) + 7.0 * boundCond(4)) * t1 +
                   15.0 * (-boundCond(3) + boundCond(0));
    nCoeffMat(2) = (0.5 * boundCond(5) - 1.5 * boundCond(2)) * t2 -
                   (6.0 * boundCond(1) + 4.0 * boundCond(4)) * t1 +
                   10.0 * (boundCond(3) - boundCond(0));
    nCoeffMat(3) = 0.5 * boundCond(2) * t2;
    nCoeffMat(4) = boundCond(1) * t1;
    nCoeffMat(5) = boundCond(0);

    double t = 1.0;
    for (int i = TrajOrder; i >= 0; i--)
    {
      CoeffMat(i) = nCoeffMat(i) / t;
      t *= dur;
    }
  }

  inline int getDim() const
  {
    return TrajDim;
  }

  inline int getOrder() const
  {
    return TrajOrder;
  }

  inline double getDuration() const
  {
    return duration;
  }

  // Get the position at time t in this piece
  inline double getPos(double t) const
  {
    // Normalize the time
    double pos(0.0);
    double tn(1.0);
    for (int i = TrajOrder; i >= 0; i--)
    {
      pos += tn * CoeffMat(i);
      tn *= t;
    }
    // The pos is not affected by normalization
    return pos;
  }

  // Get the velocity at time t in this piece
  inline double getVel(double t) const
  {
    // Normalize the time
    double vel(0.0);
    double tn(1.0);
    int n = 1;
    for (int i = TrajOrder - 1; i >= 0; i--)
    {
      vel += n * tn * CoeffMat(i);
      tn *= t;
      n++;
    }
    // Recover the actual vel
    return vel;
  }

  // Get the acceleration at time t in this piece
  inline double getAcc(double t) const
  {
    // Normalize the time
    double acc(0.0);
    double tn(1.0);
    int m = 1;
    int n = 2;
    for (int i = TrajOrder - 2; i >= 0; i--)
    {
      acc += m * n * tn * CoeffMat(i);
      tn *= t;
      m++;
      n++;
    }
    // Recover the actual acc
    return acc;
  }

  // Get the coefficient matrix of the piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline YawCoefficientMat getCoeffMat(bool normalized = false) const
  {
    YawCoefficientMat posCoeffsMat;
    double t = 1;
    for (int i = TrajOrder; i >= 0; i--)
    {
      posCoeffsMat(i) = CoeffMat(i) * t;
      t *= normalized ? duration : 1.0;
    }
    return posCoeffsMat;
  }

  // Get the polynomial coefficients of velocity of this piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline YawDotCoefficientMat getVelCoeffMat(bool normalized = false) const
  {
    YawDotCoefficientMat velCoeffMat;
    int n = 1;
    double t = 1.0;
    t *= normalized ? duration : 1.0;
    for (int i = TrajOrder - 1; i >= 0; i--)
    {
      velCoeffMat(i) = n * CoeffMat(i) * t;
      n++;
      t *= normalized ? duration : 1.0;
    }
    return velCoeffMat;
  }

  // Get the polynomial coefficients of acceleration of this piece
  // Default arg chooses the natural coefficients
  // If normalized version is needed, set the arg true
  inline YawDdotCoefficientMat getAccCoeffMat(bool normalized = false) const
  {
    YawDdotCoefficientMat accCoeffMat;
    int n = 2;
    int m = 1;
    double t = 1.0;
    t *= normalized ? duration * duration : 1.0;
    for (int i = TrajOrder - 2; i >= 0; i--)
    {
      accCoeffMat(i) = n * m * CoeffMat(i) * t;
      n++;
      m++;
      t *= normalized ? duration : 1.0;
    }
    return accCoeffMat;
  }
};

// A whole trajectory which contains multiple pieces
class Trajectory
{
private:
  typedef std::vector<Piece> Pieces;
  Pieces pieces;

public:
  Trajectory() = default;

  // Constructor from durations and coefficient matrices
  Trajectory(const std::vector<double> &durs,
             const std::vector<CoefficientMat> &coeffMats)
  {
    int N = std::min(durs.size(), coeffMats.size());
    pieces.reserve(N);
    for (int i = 0; i < N; i++)
    {
      pieces.emplace_back(durs[i], coeffMats[i]);
    }
  }

  // Constructor from pieces
  Trajectory(const Pieces &pieces) : pieces(pieces) {}

  inline int getPieceNum() const
  {
    return pieces.size();
  }

  // Get durations vector of all pieces
  inline Eigen::VectorXd getDurations() const
  {
    int N = getPieceNum();
    Eigen::VectorXd durations(N);
    for (int i = 0; i < N; i++)
    {
      durations(i) = pieces[i].getDuration();
    }
    return durations;
  }

  // Get total duration of the trajectory
  inline double getTotalDuration() const
  {
    double totalDuration = 0.0;
    for (int i = 0; i < getPieceNum(); i++)
    {
      totalDuration += pieces[i].getDuration();
    }
    return totalDuration;
  }

  // Reload the operator[] to reach the i-th piece
  inline const Piece &operator[](int i) const
  {
    return pieces[i];
  }

  inline Piece &operator[](int i)
  {
    return pieces[i];
  }

  inline void clear(void)
  {
    pieces.clear();
  }

  inline Pieces::const_iterator begin() const
  {
    return pieces.begin();
  }

  inline Pieces::const_iterator end() const
  {
    return pieces.end();
  }

  inline void reserve(const int &n)
  {
    pieces.reserve(n);
    return;
  }

  // Put another piece at the tail of this trajectory
  inline void emplace_back(const Piece &piece)
  {
    pieces.emplace_back(piece);
    return;
  }

  // Two corresponding constructors of Piece both are supported here
  template <typename ArgTypeL, typename ArgTypeR>
  inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR)
  {
    pieces.emplace_back(argL, argR);
    return;
  }

  // Append another Trajectory at the tail of this trajectory
  inline void append(const Trajectory &traj)
  {
    pieces.insert(pieces.end(), traj.begin(), traj.end());
    return;
  }

  // Append another Piece at the tail of this trajectory
  inline void append(const Piece &piece)
  {
    pieces.push_back(piece);
    return;
  }

  // Find the piece at which the time t is located
  // The index is returned and the offset in t is removed
  inline int locatePieceIdx(double &t) const
  {
    int idx;
    double dur;
    const int piece_num = getPieceNum();
    for (idx = 0;
         idx < piece_num &&
         t > (dur = pieces[idx].getDuration());
         idx++)
    {
      t -= dur;
    }
    if (idx == piece_num)
    {
      idx--;
      t += pieces[idx].getDuration();
    }
    return idx;
  }

  // Get the position at time t of the trajectory
  inline Eigen::Vector3d getPos(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getPos(t);
  }

  // Get the velocity at time t of the trajectory
  inline Eigen::Vector3d getVel(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getVel(t);
  }

  // Get the acceleration at time t of the trajectory
  inline Eigen::Vector3d getAcc(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getAcc(t);
  }

  // Get the jerk at time t of the trajectory
  inline Eigen::Vector3d getJerk(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getJerk(t);
  }

  // Get the snap at time t of the trajectory
  inline Eigen::Vector3d getSnap(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getSnap(t);
  }

  // Get the position at the juncIdx-th waypoint
  inline Eigen::Vector3d getJuncPos(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getPos(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
    }
  }

  // Get the velocity at the juncIdx-th waypoint
  inline Eigen::Vector3d getJuncVel(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getVel(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
    }
  }

  // Get the acceleration at the juncIdx-th waypoint
  inline Eigen::Vector3d getJuncAcc(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getAcc(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
    }
  }

  // Get the max velocity rate of the trajectory
  inline double getMaxVelRate() const
  {
    double maxVelRate = -INFINITY;
    double tempNorm;
    for (int i = 0; i < getPieceNum(); i++)
    {
      tempNorm = pieces[i].getMaxVelRate();
      maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
    }
    return maxVelRate;
  }

  // Get the max acceleration rate of the trajectory
  inline double getMaxAccRate() const
  {
    double maxAccRate = -INFINITY;
    double tempNorm;
    for (int i = 0; i < getPieceNum(); i++)
    {
      tempNorm = pieces[i].getMaxAccRate();
      maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
    }
    return maxAccRate;
  }

  // Get the max jerk rate of the trajectory
  inline double getMaxJerkRate() const
  {
    double maxJerkRate = -INFINITY;
    double tempNorm;
    for (int i = 0; i < getPieceNum(); i++)
    {
      tempNorm = pieces[i].getMaxJerkRate();
      maxJerkRate = maxJerkRate < tempNorm ? tempNorm : maxJerkRate;
    }
    return maxJerkRate;
  }

  // Check whether the velocity rate of this trajectory exceeds the threshold
  inline bool checkMaxVelRate(double maxVelRate) const
  {
    bool feasible = true;
    for (int i = 0; i < getPieceNum() && feasible; i++)
    {
      feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
    }
    return feasible;
  }

  // Check whether the acceleration rate of this trajectory exceeds the threshold
  inline bool checkMaxAccRate(double maxAccRate) const
  {
    bool feasible = true;
    for (int i = 0; i < getPieceNum() && feasible; i++)
    {
      feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
    }
    return feasible;
  }

  // Check whether the jerk rate of this trajectory exceeds the threshold
  inline bool checkMaxJerkRate(double maxJerkRate) const
  {
    bool feasible = true;
    for (int i = 0; i < getPieceNum() && feasible; i++)
    {
      feasible = feasible && pieces[i].checkMaxJerkRate(maxJerkRate);
    }
    return feasible;
  }

  // Scale the Trajectory(t) to Trajectory(k*t)
  inline void scaleTime(double k)
  {
    for (int i = 0; i < getPieceNum(); i++)
    {
      pieces[i].scaleTime(k);
    }
  }

  inline void sampleWholeTrajectory(std::vector<StatePVA> *vis_x) const
  {
    int n = getPieceNum();
    for (int i = 0; i < n; ++i)
    {
      pieces[i].sampleOneSeg(vis_x);
    }
  }

  inline void getWpts(std::vector<StatePVA> *wpts) const
  {
    Eigen::Vector3d pos, vel, acc;
    StatePVA x;
    pos = pieces[0].getPos(0);
    vel = pieces[0].getVel(0);
    acc = pieces[0].getAcc(0);
    x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
    wpts->push_back(x);

    int n = getPieceNum();
    for (int i = 0; i < n; ++i)
    {
      double t = pieces[i].getDuration();
      pos = pieces[i].getPos(t);
      vel = pieces[i].getVel(t);
      acc = pieces[i].getAcc(t);
      x << pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), acc(0), acc(1), acc(2);
      wpts->push_back(x);
    }
  }

  inline double getTotalDistance(const double &dt) const
  {
    double dist(0.0);
    Eigen::Vector3d p0 = getPos(0.0), p;
    for (double t = dt; t <= getTotalDuration(); t += dt)
    {
      p = getPos(t);
      dist += (p - p0).norm();
      p0 = p;
    }
    return dist;
  }

  inline double getDistance(const double &dt, const double &s, const double &e) const
  {
    double dist(0.0);
    Eigen::Vector3d p0 = getPos(s), p;
    for (double t = s + dt; t <= e; t += dt)
    {
      p = getPos(t);
      dist += (p - p0).norm();
      p0 = p;
    }
    return dist;
  }

  inline const Piece &getPiece(int i) const
  {
    return pieces[i];
  }

  inline double project_pt(const Eigen::Vector3d &pt,
                           int &ii, double &tt, Eigen::Vector3d &pro_pt)
  {
    double min_dist = -1;
    for (int i = 0; i < getPieceNum(); ++i)
    {
      auto piece = pieces[i];
      double t = 0;
      double dist = piece.project_pt(pt, t, pro_pt);
      if (dist < 0)
      {
        continue;
      }
      if (min_dist < 0 || dist < min_dist)
      {
        min_dist = dist;
        ii = i;
        tt = t;
      }
    }
    return min_dist;
  }

  inline double find_near_pt(const Eigen::Vector3d &pt, int &ii, int &ee, int &i_best, double &t_best)
  {
    double min_dist = -1;
    Eigen::Vector3d pro_pt;
    for (int i = ii; i < ee; ++i)
    {
      auto piece = pieces[i];
      double t = 0;
      double dist = piece.project_pt(pt, t, pro_pt);
      if (dist < 0)
      {
        continue;
      }
      if (min_dist < 0 || dist < min_dist)
      {
        min_dist = dist;
        i_best = i;
        t_best = t;
      }
    }
    return min_dist;
  }
};

// A whole trajectory which contains multiple pieces
class Trajectory1D
{
private:
  typedef std::vector<Piece1D> Pieces1D;
  Pieces1D pieces;

public:
  Trajectory1D() = default;

  // Constructor from durations and coefficient matrices
  Trajectory1D(const std::vector<double> &durs,
               const std::vector<YawCoefficientMat> &coeffMats)
  {
    int N = std::min(durs.size(), coeffMats.size());
    pieces.reserve(N);
    for (int i = 0; i < N; i++)
    {
      pieces.emplace_back(durs[i], coeffMats[i]);
    }
  }

  // Constructor from pieces
  Trajectory1D(const Pieces1D &pieces) : pieces(pieces) {}

  inline int getPieceNum() const
  {
    return pieces.size();
  }

  // Get durations vector of all pieces
  inline Eigen::VectorXd getDurations() const
  {
    int N = getPieceNum();
    Eigen::VectorXd durations(N);
    for (int i = 0; i < N; i++)
    {
      durations(i) = pieces[i].getDuration();
    }
    return durations;
  }

  // Get total duration of the trajectory
  inline double getTotalDuration() const
  {
    double totalDuration = 0.0;
    for (int i = 0; i < getPieceNum(); i++)
    {
      totalDuration += pieces[i].getDuration();
    }
    return totalDuration;
  }

  // Reload the operator[] to reach the i-th piece
  inline const Piece1D &operator[](int i) const
  {
    return pieces[i];
  }

  inline Piece1D &operator[](int i)
  {
    return pieces[i];
  }

  inline void clear(void)
  {
    pieces.clear();
  }

  inline Pieces1D::const_iterator begin() const
  {
    return pieces.begin();
  }

  inline Pieces1D::const_iterator end() const
  {
    return pieces.end();
  }

  inline void reserve(const int &n)
  {
    pieces.reserve(n);
    return;
  }

  // Put another piece at the tail of this trajectory
  inline void emplace_back(const Piece1D &piece)
  {
    pieces.emplace_back(piece);
    return;
  }

  // Two corresponding constructors of Piece both are supported here
  template <typename ArgTypeL, typename ArgTypeR>
  inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR)
  {
    pieces.emplace_back(argL, argR);
    return;
  }

  // Append another Trajectory at the tail of this trajectory
  inline void append(const Trajectory1D &traj)
  {
    pieces.insert(pieces.end(), traj.begin(), traj.end());
    return;
  }

  // Append another Piece at the tail of this trajectory
  inline void append(const Piece1D &piece)
  {
    pieces.push_back(piece);
    return;
  }

  // Find the piece at which the time t is located
  // The index is returned and the offset in t is removed
  inline int locatePieceIdx(double &t) const
  {
    int idx;
    double dur;
    const int piece_num = getPieceNum();
    for (idx = 0;
         idx < piece_num &&
         t > (dur = pieces[idx].getDuration());
         idx++)
    {
      t -= dur;
    }
    if (idx == piece_num)
    {
      idx--;
      t += pieces[idx].getDuration();
    }
    return idx;
  }

  // Get the position at time t of the trajectory
  inline double getPos(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getPos(t);
  }

  // Get the velocity at time t of the trajectory
  inline double getVel(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getVel(t);
  }

  // Get the acceleration at time t of the trajectory
  inline double getAcc(double t) const
  {
    int pieceIdx = locatePieceIdx(t);
    return pieces[pieceIdx].getAcc(t);
  }

  // Get the position at the juncIdx-th waypoint
  inline double getJuncPos(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getPos(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
    }
  }

  // Get the velocity at the juncIdx-th waypoint
  inline double getJuncVel(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getVel(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
    }
  }

  // Get the acceleration at the juncIdx-th waypoint
  inline double getJuncAcc(int juncIdx) const
  {
    if (juncIdx != getPieceNum())
    {
      return pieces[juncIdx].getAcc(0.0);
    }
    else
    {
      return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
    }
  }

  inline const Piece1D &getPiece(int i) const
  {
    return pieces[i];
  }
};