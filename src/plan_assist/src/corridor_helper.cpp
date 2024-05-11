#include "corridor_helper/corridor_helper.h"
using namespace std;

namespace corridor_opt
{

  bool CorridorHelper::isInhpolys(const Eigen::Vector3d pos, const Eigen::MatrixXd hpoly, double &min_d)
  {
    double minimum_d = abs(INFINITY);
    for (int i = 0; i < hpoly.cols(); i++)
    {
      Eigen::Vector3d p, n;
      p = hpoly.col(i).tail(3);
      n = hpoly.col(i).head(3);
      double d = (p - pos).dot(n);
      if (d < 1e-4)
        return false;
      if (d < minimum_d)
        minimum_d = d;
    }
    min_d = minimum_d;
    return true;
  }

  bool CorridorHelper::isInhpolysConversative(const Eigen::Vector3d pos, const Eigen::MatrixXd hpoly, double &min_d)
  {
    double minimum_d = abs(INFINITY);
    for (int i = 0; i < hpoly.cols(); i++)
    {
      Eigen::Vector3d p, n;
      p = hpoly.col(i).tail(3);
      n = hpoly.col(i).head(3);
      double d = (p - pos).dot(n);
      if (d < 0.2)
        return false;
      if (d < minimum_d)
        minimum_d = d;
    }
    min_d = minimum_d;
    return true;
  }

  int CorridorHelper::getBestHpoly(const Eigen::Vector3d pos, const int huerist_idx)
  {
    double d_best = -1;
    int idx_best = -1;
    double d;
    for (int i = huerist_idx; i < num_corridor; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        if (d > d_best)
        {
          d_best = d;
          idx_best = i;
        }
      }
      else
      {
        break;
      }
    }

    return idx_best;
  }

  int CorridorHelper::getBestHpoly(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num)
  {
    double d_best = -1;
    int idx_best = -1;
    double d;
    int index_forward = huerist_idx + expected_forward_num > num_corridor ? num_corridor : huerist_idx + expected_forward_num;
    for (int i = huerist_idx; i < index_forward; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        if (d > d_best)
        {
          d_best = d;
          idx_best = i;
        }
      }
      else
      {
        break;
      }
    }

    return idx_best;
  }

  double CorridorHelper::getBestDistInHpoly(const Eigen::Vector3d pos, const int huerist_idx)
  {
    double d_best = -1;
    double d;
    for (int i = huerist_idx; i < num_corridor; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        if (d > d_best)
        {
          d_best = d;
        }
      }
      else
      {
        break;
      }
    }

    return d_best;
  }

  double CorridorHelper::getBestDistInHpoly(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num)
  {
    double d_best = -1;
    double d;
    int index_forward = huerist_idx + expected_forward_num > num_corridor ? num_corridor : huerist_idx + expected_forward_num;
    for (int i = huerist_idx; i < index_forward; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        if (d > d_best)
        {
          d_best = d;
        }
      }
      else
      {
        break;
      }
    }
    return d_best;
  }

  int CorridorHelper::getfarthestHpolyForward(const Eigen::Vector3d pos, const int huerist_idx, const int expected_forward_num)
  {
    int idx = -1;
    double d;
    int index_forward = huerist_idx + expected_forward_num > num_corridor ? num_corridor : huerist_idx + expected_forward_num;
    for (int i = huerist_idx; i < index_forward; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        idx = i;
      }
    }

    return idx;
  }

  int CorridorHelper::getCurrentHpolyForward(const Eigen::Vector3d pos, const int huerist_idx)
  {
    int idx = -1;
    double d;
    for (int i = huerist_idx; i < num_corridor; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        idx = i;
        break;
      }
    }
    return idx;
  }

  int CorridorHelper::getCurrentHpolyBackward(const Eigen::Vector3d pos, const int huerist_idx)
  {
    int idx = -1;
    double d;
    for (int i = huerist_idx; i >= 0; i--)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        idx = i;
        break;
      }
    }
    return idx;
  }

  // only check sequence two is not enough
  int CorridorHelper::SafeCheckHpolyForward(const Eigen::Vector3d pos, const int huerist_idx)
  {
    int idx = -1;
    double d;
    int forward_id = num_corridor > huerist_idx + 10 ? huerist_idx + 10 : num_corridor;
    for (int i = huerist_idx; i < forward_id; i++)
    {
      if (isInhpolysConversative(pos, corridors_[i], d))
      {
        idx = i;
        break;
      }
    }
    return idx;
  }

  int CorridorHelper::SafeCheckHpolyBackward(const Eigen::Vector3d pos, const int huerist_idx)
  {
    int idx = -1;
    double d;
    int backward_id = huerist_idx - 3 >= 0 ? huerist_idx - 3 : 0;
    for (int i = huerist_idx; i >= backward_id; i--)
    {
      if (isInhpolysConversative(pos, corridors_[i], d))
      {
        idx = i;
        break;
      }
    }
    return idx;
  }

  int CorridorHelper::SafeCheckHpolyRange(const Eigen::Vector3d pos, const int huerist_idx, int lower, int upper)
  {
    int idx = -1;
    double d;
    int forward_id = num_corridor > huerist_idx + upper ? huerist_idx + upper : num_corridor;
    int backward_id = huerist_idx - lower >= 0 ? huerist_idx - lower : 0;
    for (int i = backward_id; i < forward_id; i++)
    {
      if (isInhpolys(pos, corridors_[i], d))
      {
        idx = i;
        break;
      }
    }
    return idx;
  }

};