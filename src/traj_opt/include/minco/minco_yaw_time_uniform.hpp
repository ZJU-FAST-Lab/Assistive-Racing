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

namespace minco
{
    /**
     * @brief
     * In this version, the middle waypoints and the segment durations are all optimizable,
     * and the segment durations are uniform.
     * min-jerk (5th degree polynomial) is provided.
     */
    namespace uniform_T
    {

        class UniformMinJerk1D
        {
        public:
            UniformMinJerk1D() = default;
            ~UniformMinJerk1D() { A.destroy(); }

        private:
            int N;
            Eigen::VectorXd headPVA, tailPVA;
            BandedSystem A;
            Eigen::VectorXd b, c;
            Eigen::VectorXd adjScaledGrad;
            Eigen::Matrix<double, 6, 1> t, tInv;


            // Eigen::MatrixXd M, M_inv;

        public : 
            inline void resetHeadCon(const Eigen::VectorXd &headState)
            {
                headPVA = headState;
            }

            inline void resetTailCon(const Eigen::VectorXd &tailState)
            {
                tailPVA = tailState;
            }

            inline void reset(const Eigen::VectorXd &headState,
                              const Eigen::VectorXd &tailState,
                              const int &pieceNum)
            {
                N = pieceNum;
                headPVA = headState;
                tailPVA = tailState;
                A.create(6 * N, 6, 6);
                b.resize(6 * N);
                c.resize(6 * N);
                adjScaledGrad.resize(6 * N);
                t(0) = 1.0;

                A(0, 0) = 1.0;
                A(1, 1) = 1.0;
                A(2, 2) = 2.0;
                for (int i = 0; i < N - 1; i++)
                {
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

            inline void generate(const Eigen::VectorXd &innerPs,
                                 const double &totalT)
            {
                t(1) = totalT / N;
                t(2) = t(1) * t(1);
                t(3) = t(2) * t(1);
                t(4) = t(2) * t(2);
                t(5) = t(4) * t(1);
                tInv = t.cwiseInverse();

                b.setZero();
                b(0) = headPVA(0);
                b(1) = headPVA(1) * t(1);
                b(2) = headPVA(2) * t(2);
                // b.row(0) = headPVA.col(0).transpose();
                // b.row(1) = headPVA.col(1).transpose() * t(1);
                // b.row(2) = headPVA.col(2).transpose() * t(2);
                for (int i = 0; i < N - 1; i++)
                {
                    b(6 * i + 5) = innerPs(i);
                }
                b(6 * N - 3) = tailPVA(0);
                b(6 * N - 2) = tailPVA(1) * t(1);
                b(6 * N - 1) = tailPVA(2) * t(2);
                // b.row(6 * N - 3) = tailPVA.col(0).transpose();
                // b.row(6 * N - 2) = tailPVA.col(1).transpose() * t(1);
                // b.row(6 * N - 1) = tailPVA.col(2).transpose() * t(2);

                A.solve(b);

                for (int i = 0; i < N; i++)
                {
                    c.block<6, 1>(6 * i, 0) =
                        b.block<6, 1>(6 * i, 0).array().colwise() * tInv.array();
                }

                return;
            }

            // inline void generate_WXX(const Eigen::Matrix3Xd &innerPs,
            //                          const double &totalT) {
            //     t(1) = totalT / N;
            //     t(2) = t(1) * t(1);
            //     t(3) = t(2) * t(1);
            //     t(4) = t(2) * t(2);
            //     t(5) = t(4) * t(1);
            //     tInv = t.cwiseInverse();

            //     b.setZero();
            //     b.row(0) = headPVA.col(0).transpose();
            //     b.row(1) = headPVA.col(1).transpose() * t(1);
            //     b.row(2) = headPVA.col(2).transpose() * t(2);
            //     for (int i = 0; i < N - 1; i++) {
            //         b.row(6 * i + 5) = innerPs.col(i).transpose();
            //     }
            //     b.row(6 * N - 3) = tailPVA.col(0).transpose();
            //     b.row(6 * N - 2) = tailPVA.col(1).transpose() * t(1);
            //     b.row(6 * N - 1) = tailPVA.col(2).transpose() * t(2);

            //     b = M_inv * b;

            //     for (int i = 0; i < N; i++) {
            //         c.block<6, 3>(6 * i, 0) =
            //             b.block<6, 3>(6 * i, 0).array().colwise() * tInv.array();
            //     }

            //     return;
            // }

            inline Trajectory1D getTraj() const
            {
                Trajectory1D traj;
                traj.reserve(N);
                for (int i = 0; i < N; i++)
                {
                    traj.emplace_back(t(1),
                                      c.block<6, 1>(6 * i, 0)
                                          .transpose()
                                          .rowwise()
                                          .reverse());
                }
                return traj;
            }

            inline double getEnergy() const
            {
                double energy(0.0);
                for (int i = 0; i < N; i++)
                {
                    energy += 36.0 * c(6 * i + 3) * c(6 * i + 3) * t(1) +
                              144.0 * c(6 * i + 4) * c(6 * i + 3) * t(2) +
                              192.0 * c(6 * i + 4) * c(6 * i + 4) * t(3) +
                              240.0 * c(6 * i + 5) * c(6 * i + 3) * t(3) +
                              720.0 * c(6 * i + 5) * c(6 * i + 4) * t(4) +
                              720.0 * c(6 * i + 5) * c(6 * i + 5) * t(5);
                }
                return energy;
            }

            inline const Eigen::VectorXd &getCoeffs(void) const
            {
                return c;
            }

            inline void getEnergyPartialGradByCoeffs(Eigen::VectorXd &gdC) const
            {
                gdC.resize(6 * N);
                for (int i = 0; i < N; i++)
                {
                    gdC(6 * i + 5) = 240.0 * c(6 * i + 3) * t(3) +
                                         720.0 * c(6 * i + 4) * t(4) +
                                         1440.0 * c(6 * i + 5) * t(5);
                    gdC(6 * i + 4) = 144.0 * c(6 * i + 3) * t(2) +
                                         384.0 * c(6 * i + 4) * t(3) +
                                         720.0 * c(6 * i + 5) * t(4);
                    gdC(6 * i + 3) = 72.0 * c(6 * i + 3) * t(1) +
                                         144.0 * c(6 * i + 4) * t(2) +
                                         240.0 * c(6 * i + 5) * t(3);
                    gdC.block<3, 1>(6 * i, 0).setZero();
                }
                return;
            }

            inline void getEnergyPartialGradByTotalTime(double &gdT) const
            {
                gdT = 0.0;
                for (int i = 0; i < N; i++)
                {
                    gdT += 36.0 * c(6 * i + 3) * c(6 * i + 3) +
                           288.0 * c(6 * i + 4) * c(6 * i + 3) * t(1) +
                           576.0 * c(6 * i + 4) * c(6 * i + 4) * t(2) +
                           720.0 * c(6 * i + 5) * c(6 * i + 3) * t(2) +
                           2880.0 * c(6 * i + 5) * c(6 * i + 4) * t(3) +
                           3600.0 * c(6 * i + 5) * c(6 * i + 5) * t(4);
                }
                gdT /= N;
                return;
            }

            inline void propogateGrad(const Eigen::VectorXd &partialGradByCoeffs,
                                      const double &partialGradByTotalTime,
                                      Eigen::VectorXd &gradByPoints,
                                      double &gradByTotalTime)
            {
                for (int i = 0; i < N; i++)
                {
                    adjScaledGrad.block<6, 1>(6 * i, 0) =
                        partialGradByCoeffs.block<6, 1>(6 * i, 0).array().colwise() * tInv.array();
                }

                A.solveAdj(adjScaledGrad);

                gradByPoints.resize(N - 1);
                for (int i = 0; i < N - 1; i++)
                {
                    gradByPoints(i) = adjScaledGrad(6 * i + 5);
                }

                gradByTotalTime = 0.0;
                gradByTotalTime += headPVA(1) * adjScaledGrad(1);
                gradByTotalTime += headPVA(2) * adjScaledGrad(2) * 2.0 * t(1);
                gradByTotalTime += tailPVA(1) * adjScaledGrad(6 * N - 2);
                gradByTotalTime += tailPVA(2) * adjScaledGrad(6 * N - 1) * 2.0 * t(1);
                Eigen::Matrix<double, 6, 1> gdtInv;
                gdtInv(0) = 0.0;
                gdtInv(1) = -1.0 * tInv(2);
                gdtInv(2) = -2.0 * tInv(3);
                gdtInv(3) = -3.0 * tInv(4);
                gdtInv(4) = -4.0 * tInv(5);
                gdtInv(5) = -5.0 * tInv(5) * tInv(1);
                // const Eigen::VectorXd gdcol = partialGradByCoeffs.cwiseProduct(b).rowwise().sum();
                const Eigen::VectorXd gdcol = partialGradByCoeffs.cwiseProduct(b);
                for (int i = 0; i < N; i++)
                {
                    gradByTotalTime += gdtInv.dot(gdcol.segment<6>(6 * i));
                }
                gradByTotalTime /= N;
                gradByTotalTime += partialGradByTotalTime;
            }

            inline void propogateGrad(const Eigen::VectorXd &partialGradByCoeffs,
                                      Eigen::VectorXd &gradByPoints)
            {
                for (int i = 0; i < N; i++)
                {
                    adjScaledGrad.block<6, 1>(6 * i, 0) =
                        partialGradByCoeffs.block<6, 1>(6 * i, 0).array().colwise() * tInv.array();
                }

                A.solveAdj(adjScaledGrad);

                gradByPoints.resize(N - 1);
                for (int i = 0; i < N - 1; i++)
                {
                    gradByPoints(i) = adjScaledGrad(6 * i + 5);
                }
            }

            inline void getPartailGradsByStartState(Eigen::VectorXd &grad)
            {
                grad = t.head<3>().asDiagonal() * adjScaledGrad.topRows<3>();
            }

            inline void getPartailGradsByEndState(Eigen::VectorXd &grad)
            {
                grad = t.head<3>().asDiagonal() * adjScaledGrad.bottomRows<3>();
            }
        };

    } // namespace uniform_T
} // namespace minco