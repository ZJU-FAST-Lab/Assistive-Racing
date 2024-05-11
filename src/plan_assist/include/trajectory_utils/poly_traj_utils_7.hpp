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

typedef Eigen::Matrix<double, 3, 8> BoundaryCond7;
typedef Eigen::Matrix<double, 3, 8> CoefficientMat7;
typedef Eigen::Matrix<double, 3, 7> VelCoefficientMat7;
typedef Eigen::Matrix<double, 3, 6> AccCoefficientMat7;
typedef Eigen::Matrix<double, 3, 5> JerCoefficientMat7;

class Piece7 {
   private:
    double duration;
    CoefficientMat7 coeffMat;

   public:
    Piece7() = default;

    Piece7(double dur, const CoefficientMat7 &cMat)
        : duration(dur), coeffMat(cMat) {}

    // Constructor from boundary condition and duration
    Piece7(BoundaryCond7 boundCond, double dur) : duration(dur) {
        double t1 = duration;
        double t2 = t1 * t1;
        double t3 = t2 * t1;

        CoefficientMat7 nCoeffMat;
        nCoeffMat.col(0) = (boundCond.col(7) / 6.0 + boundCond.col(3) / 6.0) * t3 +
                            (-2.0 * boundCond.col(6) + 2.0 * boundCond.col(2)) * t2 +
                            (10.0 * boundCond.col(5) + 10.0 * boundCond.col(1)) * t1 +
                            (-20.0 * boundCond.col(4) + 20.0 * boundCond.col(0));
        nCoeffMat.col(1) = (-0.5 * boundCond.col(7) - boundCond.col(3) / 1.5) * t3 +
                            (6.5 * boundCond.col(6) - 7.5 * boundCond.col(2)) * t2 +
                            (-34.0 * boundCond.col(5) - 36.0 * boundCond.col(1)) * t1 +
                            (70.0 * boundCond.col(4) - 70.0 * boundCond.col(0));
        nCoeffMat.col(2) = (0.5 * boundCond.col(7) + boundCond.col(3)) * t3 +
                            (-7.0 * boundCond.col(6) + 10.0 * boundCond.col(2)) * t2 +
                            (39.0 * boundCond.col(5) + 45.0 * boundCond.col(1)) * t1 +
                            (-84.0 * boundCond.col(4) + 84.0 * boundCond.col(0));
        nCoeffMat.col(3) = (-boundCond.col(7) / 6.0 - boundCond.col(3) / 1.5) * t3 +
                            (2.5 * boundCond.col(6) - 5.0 * boundCond.col(2)) * t2 +
                            (-15.0 * boundCond.col(5) - 20.0 * boundCond.col(1)) * t1 +
                            (35.0 * boundCond.col(4) - 35.0 * boundCond.col(0));
        nCoeffMat.col(4) = boundCond.col(3) * t3 / 6.0;
        nCoeffMat.col(5) = boundCond.col(2) * t2 / 2.0;
        nCoeffMat.col(6) = boundCond.col(1) * t1;
        nCoeffMat.col(7) = boundCond.col(0);
        
        double t = 1.0;
        for (int i = 7; i >= 0; i--) {
            coeffMat.col(i) = nCoeffMat.col(i) / t;
            t *= dur;
        }
    }

    Piece7(BoundaryCond7 boundCond, double dur, bool prim) : duration(dur) {
        if(!prim) std::cout<<"wrong!"<<std::endl;
        double t1 = duration;
        double t2 = t1 * t1;
        double t3 = t2 * t1;
        double t4 = t2 * t2;

        CoefficientMat7 nCoeffMat;
        nCoeffMat.col(0) = boundCond.col(7) / 42.0 * t4 - 3.0 * boundCond.col(6) / 14.0 * t3 + 6.0 * boundCond.col(5) / 7.0 * t2 -
                           10.0 * boundCond.col(4) / 7.0 * t1 + boundCond.col(3) / 14.0 * t3 + 4.0 * boundCond.col(2) / 7.0 * t2 +
                           10.0 * boundCond.col(1) / 7.0 * t1;
        nCoeffMat.col(1) = -boundCond.col(7) / 12.0 * t4 + 5.0 * boundCond.col(6) / 6.0 * t3 - 7.0 * boundCond.col(5) / 2.0 * t2 +
                           6.0 * boundCond.col(4) * t1 - boundCond.col(3) / 3.0 * t3 - 5.0 * boundCond.col(2) / 2.0 * t2 -
                           6.0 * boundCond.col(1) * t1;
        nCoeffMat.col(2) = boundCond.col(7) / 10.0 * t4 - 11.0 * boundCond.col(6) / 10.0 * t3 + 5.0 * boundCond.col(5) * t2 -
                           9.0 * boundCond.col(4) * t1 + 3.0 * boundCond.col(3) / 5.0 * t3 + 4.0 * boundCond.col(2) * t2 +
                           9.0 * boundCond.col(1) * t1;
        nCoeffMat.col(3) = -boundCond.col(7) / 24.0 * t4 + boundCond.col(6) / 2.0 *t3 - 5.0 * boundCond.col(5) / 2.0 * t2 + 
                           5.0 * boundCond.col(4) * t1 - boundCond.col(3) / 2.0 * t3 -  5.0 * boundCond.col(2) / 2.0 * t2 - 
                           5.0 * boundCond.col(1) * t1;
        nCoeffMat.col(4) = boundCond.col(3) * t3 / 6.0;
        nCoeffMat.col(5) = boundCond.col(2) * t2 / 2.0;
        nCoeffMat.col(6) = boundCond.col(1) * t1;
        nCoeffMat.col(7) = boundCond.col(0);
        
        double t = 1.0;
        for (int i = 7; i >= 0; i--) {
            coeffMat.col(i) = nCoeffMat.col(i) / t;
            t *= dur;
        }
    }

    inline int getDim() const {
        return 3;
    }

    inline int getOrder() const {
        return 7;
    }

    inline double getDuration() const {
        return duration;
    }

    inline void setDuration(double duration){
        this->duration = duration;
    }

    inline const CoefficientMat7 &getCoeffMat() const {
        return coeffMat;
    }

    inline Eigen::Vector3d getPos(const double &t) const {
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = 7; i >= 0; i--) {
            pos += tn * coeffMat.col(i);
            tn *= t;
        }
        return pos;
    }

    inline Eigen::Vector3d getVel(const double &t) const {
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = 6; i >= 0; i--) {
            vel += n * tn * coeffMat.col(i);
            tn *= t;
            n++;
        }
        return vel;
    }

    inline Eigen::Vector3d getAcc(const double &t) const {
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = 5; i >= 0; i--) {
            acc += m * n * tn * coeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        return acc;
    }

    inline Eigen::Vector3d getJer(const double &t) const {
        Eigen::Vector3d jer(0.0, 0.0, 0.0);
        double tn = 1.0;
        int l = 1;
        int m = 2;
        int n = 3;
        for (int i = 4; i >= 0; i--) {
            jer += l * m * n * tn * coeffMat.col(i);
            tn *= t;
            l++;
            m++;
            n++;
        }
        return jer;
    }

    inline Eigen::Vector3d getSna(const double &t) const {
        Eigen::Vector3d sna(0.0, 0.0, 0.0);
        double tn = 1.0;
        int l = 1;
        int m = 2;
        int n = 3;
        int o = 4;
        for (int i = 3; i >= 0; i--) {
            sna += l * m * n * o * tn * coeffMat.col(i);
            tn *= t;
            l++;
            m++;
            n++;
            o++;
        }
        return sna;
    }

    inline CoefficientMat7 normalizePosCoeffMat() const {
        CoefficientMat7 nPosCoeffsMat;
        double t = 1.0;
        for (int i = 7; i >= 0; i--) {
            nPosCoeffsMat.col(i) = coeffMat.col(i) * t;
            t *= duration;
        }
        return nPosCoeffsMat;
    }

    inline VelCoefficientMat7 normalizeVelCoeffMat() const {
        VelCoefficientMat7 nVelCoeffMat;
        int n = 1;
        double t = duration;
        for (int i = 6; i >= 0; i--) {
            nVelCoeffMat.col(i) = n * coeffMat.col(i) * t;
            t *= duration;
            n++;
        }
        return nVelCoeffMat;
    }

    inline AccCoefficientMat7 normalizeAccCoeffMat() const {
        AccCoefficientMat7 nAccCoeffMat;
        int n = 2;
        int m = 1;
        double t = duration * duration;
        for (int i = 5; i >= 0; i--) {
            nAccCoeffMat.col(i) = n * m * coeffMat.col(i) * t;
            n++;
            m++;
            t *= duration;
        }
        return nAccCoeffMat;
    }

    inline JerCoefficientMat7 normalizeJerCoeffMat() const {
        JerCoefficientMat7 nJerCoeffMat;
        int n = 3;
        int m = 2;
        int l = 1;
        double t = duration * duration * duration;
        for (int i = 4; i >= 0; i--) {
            nJerCoeffMat.col(i) = n * m * l * coeffMat.col(i) * t;
            n++;
            m++;
            l++;
            t *= duration;
        }
        return nJerCoeffMat;
    }

    inline double getMaxVelRate() const {
        VelCoefficientMat7 nVelCoeffMat = normalizeVelCoeffMat();
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return getVel(0.0).norm();
        } else {
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxVelRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    tempNormSqr = getVel((*it) * duration).squaredNorm();
                    maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                }
            }
            return sqrt(maxVelRateSqr);
        }
    }

    inline double getMaxAccRate() const {
        AccCoefficientMat7 nAccCoeffMat = normalizeAccCoeffMat();
        Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                RootFinder::polySqr(nAccCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return getAcc(0.0).norm();
        } else {
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
                r = 0.5 * (r + 1.0);
            }
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxAccRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    tempNormSqr = getAcc((*it) * duration).squaredNorm();
                    maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                }
            }
            return sqrt(maxAccRateSqr);
        }
    }

    inline double getMaxJerRate() const {
        // Compute normalized squared jerk norm polynomial coefficient matrix
        Eigen::MatrixXd nJerkCoeffMat = normalizeJerCoeffMat();
        Eigen::VectorXd coeff = RootFinder::polySqr(nJerkCoeffMat.row(0)) +
                                RootFinder::polySqr(nJerkCoeffMat.row(1)) +
                                RootFinder::polySqr(nJerkCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++) {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON) {
            return 0.0;
        } else {
            // Search an open interval whose boundaries are not zeros
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON) {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON) {
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
                 it++) {
                if (0.0 <= *it && 1.0 >= *it) {
                    // Recover the actual time then get the acc squared norm
                    tempNormSqr = getJer((*it) * duration).squaredNorm();
                    maxJerkRateSqr = maxJerkRateSqr < tempNormSqr ? tempNormSqr : maxJerkRateSqr;
                }
            }
            return sqrt(maxJerkRateSqr);
        }
    }

    inline bool checkMaxVelRate(const double &maxVelRate) const {
        double sqrMaxVelRate = maxVelRate * maxVelRate;
        if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
            getVel(duration).squaredNorm() >= sqrMaxVelRate) {
            return false;
        } else {
            VelCoefficientMat7 nVelCoeffMat = normalizeVelCoeffMat();
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(2));
            double t2 = duration * duration;
            coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    inline bool checkMaxAccRate(const double &maxAccRate) const {
        double sqrMaxAccRate = maxAccRate * maxAccRate;
        if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
            getAcc(duration).squaredNorm() >= sqrMaxAccRate) {
            return false;
        } else {
            AccCoefficientMat7 nAccCoeffMat = normalizeAccCoeffMat();
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(2));
            double t2 = duration * duration;
            double t4 = t2 * t2;
            coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    inline bool checkMaxJerRate(double maxJerkRate) const {
        double sqrMaxJerkRate = maxJerkRate * maxJerkRate;
        if (getJer(0.0).squaredNorm() >= sqrMaxJerkRate ||
            getJer(duration).squaredNorm() >= sqrMaxJerkRate) {
            return false;
        } else {
            Eigen::MatrixXd nJerkCoeffMat = normalizeJerCoeffMat();
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
};

// A whole trajectory which contains multiple pieces
class Trajectory7 {
   private:
    typedef std::vector<Piece7> Pieces;
    Pieces pieces;

   public:
    Trajectory7() = default;

    // Constructor from durations and coefficient matrices
    Trajectory7(const std::vector<double> &durs,
                const std::vector<CoefficientMat7> &coeffMats) {
        int N = std::min(durs.size(), coeffMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; i++) {
            pieces.emplace_back(durs[i], coeffMats[i]);
        }
    }

    // Constructor from durations and coefficient matrices
    Trajectory7(const Eigen::VectorXd &durs,
                const std::vector<CoefficientMat7> &coeffMats) {
        size_t n = durs.rows();
        int N = std::min(n, coeffMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; i++) {
            pieces.emplace_back(durs[i], coeffMats[i]);
        }
    }

    // Constructor from pieces
    Trajectory7(const Pieces &pieces) : pieces(pieces) {}

    inline int getPieceNum() const {
        return pieces.size();
    }

    // Get durations vector of all pieces
    inline Eigen::VectorXd getDurations() const {
        int N = getPieceNum();
        Eigen::VectorXd durations(N);
        for (int i = 0; i < N; i++) {
            durations(i) = pieces[i].getDuration();
        }
        return durations;
    }

    // Get total duration of the trajectory
    inline double getTotalDuration() const {
        double totalDuration = 0.0;
        for (int i = 0; i < getPieceNum(); i++) {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    // Reload the operator[] to reach the i-th piece
    inline const Piece7 &operator[](int i) const {
        return pieces[i];
    }

    inline Piece7 &operator[](int i) {
        return pieces[i];
    }

    inline void clear(void) {
        pieces.clear();
    }

    inline Pieces::const_iterator begin() const {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const {
        return pieces.end();
    }

    inline void reserve(const int &n) {
        pieces.reserve(n);
        return;
    }

    // Put another piece at the tail of this trajectory
    inline void emplace_back(const Piece7 &piece) {
        pieces.emplace_back(piece);
        return;
    }

    // Two corresponding constructors of Piece7 both are supported here
    template <typename ArgTypeL, typename ArgTypeR>
    inline void emplace_back(const ArgTypeL &argL, const ArgTypeR &argR) {
        pieces.emplace_back(argL, argR);
        return;
    }

    // Append another Trajectory7 at the tail of this trajectory
    inline void append(const Trajectory7 &traj) {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    // Find the piece at which the time t is located
    // The index is returned and the offset in t is removed
    inline int locatePieceIdx(double &t) const {
        int idx;
        double dur;
        const int piece_num = getPieceNum();
        for (idx = 0;
             idx < piece_num &&
             t > (dur = pieces[idx].getDuration());
             idx++) {
            t -= dur;
        }
        if (idx == piece_num) {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    // Get the position at time t of the trajectory
    inline Eigen::Vector3d getPos(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    // Get the velocity at time t of the trajectory
    inline Eigen::Vector3d getVel(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    // Get the acceleration at time t of the trajectory
    inline Eigen::Vector3d getAcc(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    // Get the jerk at time t of the trajectory
    inline Eigen::Vector3d getJerk(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getJer(t);
    }

    // Get the smap at time t of the trajectory
    inline Eigen::Vector3d getSnap(double t) const {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getSna(t);
    }

    // Get the position at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncPos(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getPos(0.0);
        } else {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the velocity at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncVel(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getVel(0.0);
        } else {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the acceleration at the juncIdx-th waypoint
    inline Eigen::Vector3d getJuncAcc(int juncIdx) const {
        if (juncIdx != getPieceNum()) {
            return pieces[juncIdx].getAcc(0.0);
        } else {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }

    // Get the max velocity rate of the trajectory
    inline double getMaxVelRate() const {
        double maxVelRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxVelRate();
            maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
        }
        return maxVelRate;
    }

    // Get the max acceleration rate of the trajectory
    inline double getMaxAccRate() const {
        double maxAccRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxAccRate();
            maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
        }
        return maxAccRate;
    }

    // Get the max jerk rate of the trajectory
    inline double getMaxJerkRate() const {
        double maxJerkRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < getPieceNum(); i++) {
            tempNorm = pieces[i].getMaxJerRate();
            maxJerkRate = maxJerkRate < tempNorm ? tempNorm : maxJerkRate;
        }
        return maxJerkRate;
    }

    // Check whether the velocity rate of this trajectory exceeds the threshold
    inline bool checkMaxVelRate(double maxVelRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
        }
        return feasible;
    }

    // Check whether the acceleration rate of this trajectory exceeds the threshold
    inline bool checkMaxAccRate(double maxAccRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
        }
        return feasible;
    }

    // Check whether the jerk rate of this trajectory exceeds the threshold
    inline bool checkMaxJerkRate(double maxJerkRate) const {
        bool feasible = true;
        for (int i = 0; i < getPieceNum() && feasible; i++) {
            feasible = feasible && pieces[i].checkMaxJerRate(maxJerkRate);
        }
        return feasible;
    }

    inline const Piece7 &getPiece(int i) const {
        return pieces[i];
    }
};