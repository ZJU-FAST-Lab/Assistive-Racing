/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)
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

#include "banded_system.hpp"
#include "trajectory_utils/poly_traj_utils.hpp"
#include "trajectory_utils/poly_traj_utils_7.hpp"

namespace minco {
/**
 * @brief 
 * In this version, the middle waypoints and the segment durations are all optimizable, 
 * and the segment durations are uniform.
 * Two kinds of control are provided: min-jerk (5th degree polynomial) and min-snap (7th).
 */
namespace uniform_T {

class UniformMinJerk {
   public:
    UniformMinJerk() = default;
    ~UniformMinJerk() { A.destroy(); }

   private:
    int N;
    Eigen::Matrix3d headPVA, tailPVA;
    BandedSystem A;
    Eigen::MatrixX3d b, c;
    Eigen::MatrixX3d adjScaledGrad;
    Eigen::Matrix<double, 6, 1> t, tInv;

    // Eigen::MatrixXd M, M_inv;

   public:
    inline void resetHeadCon(const Eigen::Matrix3d &headState) {
        headPVA = headState;
    }

    inline void resetTailCon(const Eigen::Matrix3d &tailState) {
        tailPVA = tailState;
    }

    inline void reset(const Eigen::Matrix3d &headState,
                      const Eigen::Matrix3d &tailState,
                      const int &pieceNum) {
        N = pieceNum;
        headPVA = headState;
        tailPVA = tailState;
        A.create(6 * N, 6, 6);
        b.resize(6 * N, 3);
        c.resize(6 * N, 3);
        adjScaledGrad.resize(6 * N, 3);
        t(0) = 1.0;

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        for (int i = 0; i < N - 1; i++) {
            A(6 * i + 3, 6 * i + 3) = 6.0;
            A(6 * i + 3, 6 * i + 4) = 24.0;
            A(6 * i + 3, 6 * i + 5) = 60.0;
            A(6 * i + 3, 6 * i + 9) = -6.0;
            A(6 * i + 4, 6 * i + 4) = 24.0;
            A(6 * i + 4, 6 * i + 5) = 120.0;
            A(6 * i + 4, 6 * i + 10) = -24.0;
            A(6 * i + 5, 6 * i) = 1.0;
            A(6 * i + 5, 6 * i + 1) = 1.0;
            A(6 * i + 5, 6 * i + 2) = 1.0;
            A(6 * i + 5, 6 * i + 3) = 1.0;
            A(6 * i + 5, 6 * i + 4) = 1.0;
            A(6 * i + 5, 6 * i + 5) = 1.0;
            A(6 * i + 6, 6 * i) = 1.0;
            A(6 * i + 6, 6 * i + 1) = 1.0;
            A(6 * i + 6, 6 * i + 2) = 1.0;
            A(6 * i + 6, 6 * i + 3) = 1.0;
            A(6 * i + 6, 6 * i + 4) = 1.0;
            A(6 * i + 6, 6 * i + 5) = 1.0;
            A(6 * i + 6, 6 * i + 6) = -1.0;
            A(6 * i + 7, 6 * i + 1) = 1.0;
            A(6 * i + 7, 6 * i + 2) = 2.0;
            A(6 * i + 7, 6 * i + 3) = 3.0;
            A(6 * i + 7, 6 * i + 4) = 4.0;
            A(6 * i + 7, 6 * i + 5) = 5.0;
            A(6 * i + 7, 6 * i + 7) = -1.0;
            A(6 * i + 8, 6 * i + 2) = 2.0;
            A(6 * i + 8, 6 * i + 3) = 6.0;
            A(6 * i + 8, 6 * i + 4) = 12.0;
            A(6 * i + 8, 6 * i + 5) = 20.0;
            A(6 * i + 8, 6 * i + 8) = -2.0;
        }
        A(6 * N - 3, 6 * N - 6) = 1.0;
        A(6 * N - 3, 6 * N - 5) = 1.0;
        A(6 * N - 3, 6 * N - 4) = 1.0;
        A(6 * N - 3, 6 * N - 3) = 1.0;
        A(6 * N - 3, 6 * N - 2) = 1.0;
        A(6 * N - 3, 6 * N - 1) = 1.0;
        A(6 * N - 2, 6 * N - 5) = 1.0;
        A(6 * N - 2, 6 * N - 4) = 2.0;
        A(6 * N - 2, 6 * N - 3) = 3.0;
        A(6 * N - 2, 6 * N - 2) = 4.0;
        A(6 * N - 2, 6 * N - 1) = 5.0;
        A(6 * N - 1, 6 * N - 4) = 2.0;
        A(6 * N - 1, 6 * N - 3) = 6.0;
        A(6 * N - 1, 6 * N - 2) = 12.0;
        A(6 * N - 1, 6 * N - 1) = 20.0;
        A.factorizeLU();

        return;
    }

    inline void generate(const Eigen::Matrix3Xd &innerPs,
                         const double &totalT) {
        t(1) = totalT / N;
        t(2) = t(1) * t(1);
        t(3) = t(2) * t(1);
        t(4) = t(2) * t(2);
        t(5) = t(4) * t(1);
        tInv = t.cwiseInverse();

        b.setZero();
        b.row(0) = headPVA.col(0).transpose();
        b.row(1) = headPVA.col(1).transpose() * t(1);
        b.row(2) = headPVA.col(2).transpose() * t(2);
        for (int i = 0; i < N - 1; i++) {
            b.row(6 * i + 5) = innerPs.col(i).transpose();
        }
        b.row(6 * N - 3) = tailPVA.col(0).transpose();
        b.row(6 * N - 2) = tailPVA.col(1).transpose() * t(1);
        b.row(6 * N - 1) = tailPVA.col(2).transpose() * t(2);

        A.solve(b);

        for (int i = 0; i < N; i++) {
            c.block<6, 3>(6 * i, 0) =
                b.block<6, 3>(6 * i, 0).array().colwise() * tInv.array();
        }

        return;
    }

    inline Trajectory getTraj() const {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++) {
            traj.emplace_back(t(1),
                              c.block<6, 3>(6 * i, 0)
                                  .transpose()
                                  .rowwise()
                                  .reverse());
        }
        return traj;
    }

    inline double getEnergy() const {
        double energy(0.0);
        for (int i = 0; i < N; i++) {
            energy += 36.0 * c.row(6 * i + 3).squaredNorm() * t(1) +
                      144.0 * c.row(6 * i + 4).dot(c.row(6 * i + 3)) * t(2) +
                      192.0 * c.row(6 * i + 4).squaredNorm() * t(3) +
                      240.0 * c.row(6 * i + 5).dot(c.row(6 * i + 3)) * t(3) +
                      720.0 * c.row(6 * i + 5).dot(c.row(6 * i + 4)) * t(4) +
                      720.0 * c.row(6 * i + 5).squaredNorm() * t(5);
        }
        return energy;
    }

    inline const Eigen::MatrixX3d &getCoeffs(void) const {
        return c;
    }

    inline void getEnergyPartialGradByCoeffs(Eigen::MatrixX3d &gdC) const {
        gdC.resize(6 * N, 3);
        for (int i = 0; i < N; i++) {
            gdC.row(6 * i + 5) = 240.0 * c.row(6 * i + 3) * t(3) +
                                 720.0 * c.row(6 * i + 4) * t(4) +
                                 1440.0 * c.row(6 * i + 5) * t(5);
            gdC.row(6 * i + 4) = 144.0 * c.row(6 * i + 3) * t(2) +
                                 384.0 * c.row(6 * i + 4) * t(3) +
                                 720.0 * c.row(6 * i + 5) * t(4);
            gdC.row(6 * i + 3) = 72.0 * c.row(6 * i + 3) * t(1) +
                                 144.0 * c.row(6 * i + 4) * t(2) +
                                 240.0 * c.row(6 * i + 5) * t(3);
            gdC.block<3, 3>(6 * i, 0).setZero();
        }
        return;
    }

    inline void getEnergyPartialGradByTotalTime(double &gdT) const {
        gdT = 0.0;
        for (int i = 0; i < N; i++) {
            gdT += 36.0 * c.row(6 * i + 3).squaredNorm() +
                   288.0 * c.row(6 * i + 4).dot(c.row(6 * i + 3)) * t(1) +
                   576.0 * c.row(6 * i + 4).squaredNorm() * t(2) +
                   720.0 * c.row(6 * i + 5).dot(c.row(6 * i + 3)) * t(2) +
                   2880.0 * c.row(6 * i + 5).dot(c.row(6 * i + 4)) * t(3) +
                   3600.0 * c.row(6 * i + 5).squaredNorm() * t(4);
        }
        gdT /= N;
        return;
    }

    inline void propogateGrad(const Eigen::MatrixX3d &partialGradByCoeffs,
                              const double &partialGradByTotalTime,
                              Eigen::Matrix3Xd &gradByPoints,
                              double &gradByTotalTime) {
        for (int i = 0; i < N; i++) {
            adjScaledGrad.block<6, 3>(6 * i, 0) =
                partialGradByCoeffs.block<6, 3>(6 * i, 0).array().colwise() * tInv.array();
        }

        A.solveAdj(adjScaledGrad);

        gradByPoints.resize(3, N - 1);
        for (int i = 0; i < N - 1; i++) {
            gradByPoints.col(i) = adjScaledGrad.row(6 * i + 5).transpose();
        }

        gradByTotalTime = 0.0;
        gradByTotalTime += headPVA.col(1).dot(adjScaledGrad.row(1));
        gradByTotalTime += headPVA.col(2).dot(adjScaledGrad.row(2)) * 2.0 * t(1);
        gradByTotalTime += tailPVA.col(1).dot(adjScaledGrad.row(6 * N - 2));
        gradByTotalTime += tailPVA.col(2).dot(adjScaledGrad.row(6 * N - 1)) * 2.0 * t(1);
        Eigen::Matrix<double, 6, 1> gdtInv;
        gdtInv(0) = 0.0;
        gdtInv(1) = -1.0 * tInv(2);
        gdtInv(2) = -2.0 * tInv(3);
        gdtInv(3) = -3.0 * tInv(4);
        gdtInv(4) = -4.0 * tInv(5);
        gdtInv(5) = -5.0 * tInv(5) * tInv(1);
        const Eigen::VectorXd gdcol = partialGradByCoeffs.cwiseProduct(b).rowwise().sum();
        for (int i = 0; i < N; i++) {
            gradByTotalTime += gdtInv.dot(gdcol.segment<6>(6 * i));
        }
        gradByTotalTime /= N;
        gradByTotalTime += partialGradByTotalTime;
    }

    inline void propogateGrad(const Eigen::MatrixX3d &partialGradByCoeffs,
                              Eigen::Matrix3Xd &gradByPoints){
        for (int i = 0; i < N; i++) {
            adjScaledGrad.block<6, 3>(6 * i, 0) =
                partialGradByCoeffs.block<6, 3>(6 * i, 0).array().colwise() * tInv.array();
        }

        A.solveAdj(adjScaledGrad);

        gradByPoints.resize(3, N - 1);
        for (int i = 0; i < N - 1; i++) {
            gradByPoints.col(i) = adjScaledGrad.row(6 * i + 5).transpose();
        }                            
    }

    inline void getPartailGradsByStartState(Eigen::MatrixXd &grad) {
        grad = t.head<3>().asDiagonal() * adjScaledGrad.topRows<3>();
    }

    inline void getPartailGradsByEndState(Eigen::MatrixXd &grad) {
        grad = t.head<3>().asDiagonal() * adjScaledGrad.bottomRows<3>();
    }
};

class UniformMinSnap {
   public:
    UniformMinSnap() = default;
    ~UniformMinSnap() { A.destroy(); }

   private:
    int N;
    Eigen::Matrix<double, 3, 4> headPVAJ;
    Eigen::Matrix<double, 3, 4> tailPVAJ;
    BandedSystem A;
    Eigen::MatrixX3d b;
    /**
     * @brief 
     * c = [x, y, z]
     * x = [c1^T, c2^T, ..., cN^T]^T
     * c1 = [c10, c11, c12, c13, c14, c15, c16, c17]^T
     */
    Eigen::MatrixX3d c;
    Eigen::MatrixX3d adjScaledGrad;
    /**
     * @brief 
     * t = [1, T, T^2, ..., T^7]
     */
    Eigen::Matrix<double, 8, 1> t, tInv;

   public:
    /**
    * @brief 
    * Reset the end state in case of optimizable tail condition. 
    * @param tailState [p, v, a, j], p/v/a/j = [x, y, z]^T
    */
    inline void resetTailCon(const Eigen::Matrix<double, 3, 4> &tailState) {
        tailPVAJ = tailState;
    }

    /**
     * @brief 
     * Resize each matrix, construct the UNIFIED M (T--->1) and LU factorize it.
     * @param headState [p, v, a, j], p/v/a/j = [x, y, z]^T
     * @param tailState [p, v, a, j], p/v/a/j = [x, y, z]^T
     * @param pieceNum number of segments
     * @note This function need to be called before generate()
     */
    inline void reset(const Eigen::Matrix<double, 3, 4> &headState,
                      const Eigen::Matrix<double, 3, 4> &tailState,
                      const int &pieceNum) {
        N = pieceNum;
        headPVAJ = headState;
        tailPVAJ = tailState;
        A.create(8 * N, 8, 8);
        b.resize(8 * N, 3);
        c.resize(8 * N, 3);
        adjScaledGrad.resize(8 * N, 3);
        t(0) = 1.0;

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        A(3, 3) = 6.0;
        for (int i = 0; i < N - 1; i++) {
            A(8 * i + 4, 8 * i + 4) = 24.0;
            A(8 * i + 4, 8 * i + 5) = 120.0;
            A(8 * i + 4, 8 * i + 6) = 360.0;
            A(8 * i + 4, 8 * i + 7) = 840.0;
            A(8 * i + 4, 8 * i + 12) = -24.0;

            A(8 * i + 5, 8 * i + 5) = 120.0;
            A(8 * i + 5, 8 * i + 6) = 720.0;
            A(8 * i + 5, 8 * i + 7) = 2520.0;
            A(8 * i + 5, 8 * i + 13) = -120.0;

            A(8 * i + 6, 8 * i + 6) = 720.0;
            A(8 * i + 6, 8 * i + 7) = 5040.0;
            A(8 * i + 6, 8 * i + 14) = -720.0;

            A(8 * i + 7, 8 * i) = 1.0;
            A(8 * i + 7, 8 * i + 1) = 1.0;
            A(8 * i + 7, 8 * i + 2) = 1.0;
            A(8 * i + 7, 8 * i + 3) = 1.0;
            A(8 * i + 7, 8 * i + 4) = 1.0;
            A(8 * i + 7, 8 * i + 5) = 1.0;
            A(8 * i + 7, 8 * i + 6) = 1.0;
            A(8 * i + 7, 8 * i + 7) = 1.0;

            A(8 * i + 8, 8 * i) = 1.0;
            A(8 * i + 8, 8 * i + 1) = 1.0;
            A(8 * i + 8, 8 * i + 2) = 1.0;
            A(8 * i + 8, 8 * i + 3) = 1.0;
            A(8 * i + 8, 8 * i + 4) = 1.0;
            A(8 * i + 8, 8 * i + 5) = 1.0;
            A(8 * i + 8, 8 * i + 6) = 1.0;
            A(8 * i + 8, 8 * i + 7) = 1.0;
            A(8 * i + 8, 8 * i + 8) = -1.0;

            A(8 * i + 9, 8 * i + 1) = 1.0;
            A(8 * i + 9, 8 * i + 2) = 2.0;
            A(8 * i + 9, 8 * i + 3) = 3.0;
            A(8 * i + 9, 8 * i + 4) = 4.0;
            A(8 * i + 9, 8 * i + 5) = 5.0;
            A(8 * i + 9, 8 * i + 6) = 6.0;
            A(8 * i + 9, 8 * i + 7) = 7.0;
            A(8 * i + 9, 8 * i + 9) = -1.0;

            A(8 * i + 10, 8 * i + 2) = 2.0;
            A(8 * i + 10, 8 * i + 3) = 6.0;
            A(8 * i + 10, 8 * i + 4) = 12.0;
            A(8 * i + 10, 8 * i + 5) = 20.0;
            A(8 * i + 10, 8 * i + 6) = 30.0;
            A(8 * i + 10, 8 * i + 7) = 42.0;
            A(8 * i + 10, 8 * i + 10) = -2.0;

            A(8 * i + 11, 8 * i + 3) = 6.0;
            A(8 * i + 11, 8 * i + 4) = 24.0;
            A(8 * i + 11, 8 * i + 5) = 60.0;
            A(8 * i + 11, 8 * i + 6) = 120.0;
            A(8 * i + 11, 8 * i + 7) = 210.0;
            A(8 * i + 11, 8 * i + 11) = -6.0;
        }
        A(8 * N - 4, 8 * N - 8) = 1.0;
        A(8 * N - 4, 8 * N - 7) = 1.0;
        A(8 * N - 4, 8 * N - 6) = 1.0;
        A(8 * N - 4, 8 * N - 5) = 1.0;
        A(8 * N - 4, 8 * N - 4) = 1.0;
        A(8 * N - 4, 8 * N - 3) = 1.0;
        A(8 * N - 4, 8 * N - 2) = 1.0;
        A(8 * N - 4, 8 * N - 1) = 1.0;

        A(8 * N - 3, 8 * N - 7) = 1.0;
        A(8 * N - 3, 8 * N - 6) = 2.0;
        A(8 * N - 3, 8 * N - 5) = 3.0;
        A(8 * N - 3, 8 * N - 4) = 4.0;
        A(8 * N - 3, 8 * N - 3) = 5.0;
        A(8 * N - 3, 8 * N - 2) = 6.0;
        A(8 * N - 3, 8 * N - 1) = 7.0;

        A(8 * N - 2, 8 * N - 6) = 2.0;
        A(8 * N - 2, 8 * N - 5) = 6.0;
        A(8 * N - 2, 8 * N - 4) = 12.0;
        A(8 * N - 2, 8 * N - 3) = 20.0;
        A(8 * N - 2, 8 * N - 2) = 30.0;
        A(8 * N - 2, 8 * N - 1) = 42.0;

        A(8 * N - 1, 8 * N - 5) = 6.0;
        A(8 * N - 1, 8 * N - 4) = 24.0;
        A(8 * N - 1, 8 * N - 3) = 60.0;
        A(8 * N - 1, 8 * N - 2) = 120.0;
        A(8 * N - 1, 8 * N - 1) = 210.0;

        A.factorizeLU();

        return;
    }

    /**
     * @brief 
     * Given middle waypoints (3D) and total traj duration, 
     * solve the minimum control optimization, 
     * and store the UN-UNIFIED coeffs in c.
     * @param innerPs middle waypoints, [P1, P2, ..., PN-1], P = [px, py, pz]^T
     * @param totalT total traj duration
     * @note After this operation, the c stores the UN-UNIFIED coeffs. 
     * When constructing b, it is scaled as p-->p, v-->vT, a-->aT^2, j-->jT^3, s-->sT^4 
     */
    inline void generate(const Eigen::Matrix3Xd &innerPs,
                         const double &totalT) {
        t(1) = totalT / N;
        t(2) = t(1) * t(1);
        t(3) = t(2) * t(1);
        t(4) = t(2) * t(2);
        t(5) = t(4) * t(1);
        t(6) = t(5) * t(1);
        t(7) = t(6) * t(1);
        tInv = t.cwiseInverse();

        b.setZero();
        b.row(0) = headPVAJ.col(0).transpose();
        b.row(1) = headPVAJ.col(1).transpose() * t(1);
        b.row(2) = headPVAJ.col(2).transpose() * t(2);
        b.row(3) = headPVAJ.col(3).transpose() * t(3);
        for (int i = 0; i < N - 1; i++) {
            b.row(8 * i + 7) = innerPs.col(i).transpose();
        }
        b.row(8 * N - 4) = tailPVAJ.col(0).transpose();
        b.row(8 * N - 3) = tailPVAJ.col(1).transpose() * t(1);
        b.row(8 * N - 2) = tailPVAJ.col(2).transpose() * t(2);
        b.row(8 * N - 1) = tailPVAJ.col(3).transpose() * t(3);

        A.solve(b);

        for (int i = 0; i < N; i++) {
            c.block<8, 3>(8 * i, 0) =
                b.block<8, 3>(8 * i, 0).array().colwise() * tInv.array();
        }

        return;
    }

    /**
     * @brief Get the Traj object
     * 
     * @return the minimum control 7th degree Trajectory
     */
    inline Trajectory7 getTraj() const {
        Trajectory7 traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++) {
            traj.emplace_back(t(1),
                              c.block<8, 3>(8 * i, 0)
                                  .transpose()
                                  .rowwise()
                                  .reverse());
        }
        return traj;
    }

    /**
     * @brief
     * Get the cost of control.
     * @return double 
     */
    inline double getEnergy() const {
        double energy = 0.0;
        for (int i = 0; i < N; i++) {
            energy += 576.0 * c.row(8 * i + 4).squaredNorm() * t(1) +
                      2880.0 * c.row(8 * i + 4).dot(c.row(8 * i + 5)) * t(2) +
                      4800.0 * c.row(8 * i + 5).squaredNorm() * t(3) +
                      5760.0 * c.row(8 * i + 4).dot(c.row(8 * i + 6)) * t(3) +
                      21600.0 * c.row(8 * i + 5).dot(c.row(8 * i + 6)) * t(4) +
                      10080.0 * c.row(8 * i + 4).dot(c.row(8 * i + 7)) * t(4) +
                      25920.0 * c.row(8 * i + 6).squaredNorm() * t(5) +
                      40320.0 * c.row(8 * i + 5).dot(c.row(8 * i + 7)) * t(5) +
                      100800.0 * c.row(8 * i + 6).dot(c.row(8 * i + 7)) * t(6) +
                      100800.0 * c.row(8 * i + 7).squaredNorm() * t(7);
        }
        return energy;
    }

    /**
     * @brief Get the UN-UNIFIED Coeffs
     * 
     * @return 
     * c = [x, y, z], 
     * x = [c1^T, c2^T, ..., cN^T]^T, 
     * c1 = [c10, c11, c12, c13, c14, c15, c16, c17]^T
     */
    inline const Eigen::MatrixX3d &getCoeffs(void) const {
        return c;
    }

    /**
     * @brief Get the Energy Partial Grad By Coeffs
     * 
     * @param gdC returned UN-UNIFIED grad, elements sequence corresponding to c. 
     * @note This function has to be called after generate(). 
     */
    inline void getEnergyPartialGradByCoeffs(Eigen::MatrixX3d &gdC) const {
        gdC.resize(8 * N, 3);
        for (int i = 0; i < N; i++) {
            gdC.row(8 * i + 7) = 10080.0 * c.row(8 * i + 4) * t(4) +
                                 40320.0 * c.row(8 * i + 5) * t(5) +
                                 100800.0 * c.row(8 * i + 6) * t(6) +
                                 201600.0 * c.row(8 * i + 7) * t(7);
            gdC.row(8 * i + 6) = 5760.0 * c.row(8 * i + 4) * t(3) +
                                 21600.0 * c.row(8 * i + 5) * t(4) +
                                 51840.0 * c.row(8 * i + 6) * t(5) +
                                 100800.0 * c.row(8 * i + 7) * t(6);
            gdC.row(8 * i + 5) = 2880.0 * c.row(8 * i + 4) * t(2) +
                                 9600.0 * c.row(8 * i + 5) * t(3) +
                                 21600.0 * c.row(8 * i + 6) * t(4) +
                                 40320.0 * c.row(8 * i + 7) * t(5);
            gdC.row(8 * i + 4) = 1152.0 * c.row(8 * i + 4) * t(1) +
                                 2880.0 * c.row(8 * i + 5) * t(2) +
                                 5760.0 * c.row(8 * i + 6) * t(3) +
                                 10080.0 * c.row(8 * i + 7) * t(4);
            gdC.block<4, 3>(8 * i, 0).setZero();
        }
        return;
    }

    /**
     * @brief Get the Energy Partial Grad By Total Time
     * 
     * @param gdT returned grad to toal time duration
     * @note This function has to be called after generate(), 
     * and maybe along with getEnergyPartialGradByCoeffs()
     */
    inline void getEnergyPartialGradByTotalTime(double &gdT) const {
        gdT = 0.0;
        for (int i = 0; i < N; i++) {
            gdT += 576.0 * c.row(8 * i + 4).squaredNorm() +
                   5760.0 * c.row(8 * i + 4).dot(c.row(8 * i + 5)) * t(1) +
                   14400.0 * c.row(8 * i + 5).squaredNorm() * t(2) +
                   17280.0 * c.row(8 * i + 4).dot(c.row(8 * i + 6)) * t(2) +
                   86400.0 * c.row(8 * i + 5).dot(c.row(8 * i + 6)) * t(3) +
                   40320.0 * c.row(8 * i + 4).dot(c.row(8 * i + 7)) * t(3) +
                   129600.0 * c.row(8 * i + 6).squaredNorm() * t(4) +
                   201600.0 * c.row(8 * i + 5).dot(c.row(8 * i + 7)) * t(4) +
                   604800.0 * c.row(8 * i + 6).dot(c.row(8 * i + 7)) * t(5) +
                   705600.0 * c.row(8 * i + 7).squaredNorm() * t(6);
        }
        gdT /= N;
        return;
    }

    inline void propogateGrad(const Eigen::MatrixX3d &partialGradByCoeffs,
                              Eigen::Matrix3Xd &gradByPoints){
        for (int i = 0; i < N; i++) {
            adjScaledGrad.block<8, 3>(8 * i, 0) =
                partialGradByCoeffs.block<8, 3>(8 * i, 0).array().colwise() * tInv.array();
        }

        A.solveAdj(adjScaledGrad);

        gradByPoints.resize(3, N - 1);
        for (int i = 0; i < N - 1; i++) {
            gradByPoints.col(i) = adjScaledGrad.row(8 * i + 7).transpose();
        }                            
    }

    /**
     * @brief Propogate the partial gradient by c and by total T to by middle wpts and by total T. 
     * 
     * @param partialGradByCoeffs 
     * the partial gradient of energy cost AND OTHER user-defined cost by c.
     * @param partialGradByTotalTime 
     * the partial gradient of energy cost AND OTHER user-defined cost by total T.
     * @param gradByPoints 
     * After propogation, the partial gradient of energy cost AND OTHER user-defined cost by middle wpts.
     * @param gradByTotalTime 
     * After propogation, the partial gradient of energy cost AND OTHER user-defined cost by total T.
     */
    inline void propogateGrad(const Eigen::MatrixX3d &partialGradByCoeffs,
                              const double &partialGradByTotalTime,
                              Eigen::Matrix3Xd &gradByPoints,
                              double &gradByTotalTime) {
        for (int i = 0; i < N; i++) {
            adjScaledGrad.block<8, 3>(8 * i, 0) =
                partialGradByCoeffs.block<8, 3>(8 * i, 0).array().colwise() * tInv.array();
        }

        // A is normalized before, so adjScaledGrad needs to be normalized before cal Ax=b
        A.solveAdj(adjScaledGrad);

        gradByPoints.resize(3, N - 1);
        for (int i = 0; i < N - 1; i++) {
            gradByPoints.col(i) = adjScaledGrad.row(8 * i + 7).transpose();
        }

        gradByTotalTime = 0.0;
        gradByTotalTime += headPVAJ.col(1).dot(adjScaledGrad.row(1));
        gradByTotalTime += headPVAJ.col(2).dot(adjScaledGrad.row(2)) * 2.0 * t(1);
        gradByTotalTime += headPVAJ.col(3).dot(adjScaledGrad.row(3)) * 3.0 * t(2);
        gradByTotalTime += tailPVAJ.col(1).dot(adjScaledGrad.row(8 * N - 3));
        gradByTotalTime += tailPVAJ.col(2).dot(adjScaledGrad.row(8 * N - 2)) * 2.0 * t(1);
        gradByTotalTime += tailPVAJ.col(3).dot(adjScaledGrad.row(8 * N - 1)) * 3.0 * t(2);
        Eigen::Matrix<double, 8, 1> gdtInv;
        gdtInv(0) = 0.0;
        gdtInv(1) = -1.0 * tInv(2);
        gdtInv(2) = -2.0 * tInv(3);
        gdtInv(3) = -3.0 * tInv(4);
        gdtInv(4) = -4.0 * tInv(5);
        gdtInv(5) = -5.0 * tInv(6);
        gdtInv(6) = -6.0 * tInv(7);
        gdtInv(7) = -7.0 * tInv(7) * tInv(1);
        const Eigen::VectorXd gdcol = partialGradByCoeffs.cwiseProduct(b).rowwise().sum();
        for (int i = 0; i < N; i++) {
            gradByTotalTime += gdtInv.dot(gdcol.segment<8>(8 * i));
        }
        gradByTotalTime /= N;
        gradByTotalTime += partialGradByTotalTime;
    }

    /**
     * @brief Get the Partial Grad by end state
     * 
     * @param grad the partial grad by end state. 
     * end state(3*4) = [p v a j]^T, 
     * p = [x, y, z]^T.
     */
    inline void getPartailGradsByEndState(Eigen::MatrixXd &grad) {
        grad = t.head<4>().asDiagonal() * adjScaledGrad.bottomRows<4>();
    }
};

}  // namespace uniform_T
}  // namespace minco