Header header
geometry_msgs/Point position
geometry_msgs/Vector3 velocity
geometry_msgs/Vector3 acceleration
geometry_msgs/Vector3 jerk
float64 yaw
float64 yaw_dot
float64[3] kx
float64[3] kv 

uint32 trajectory_id

uint8 TRAJECTORY_STATUS_EMPTY = 0
uint8 TRAJECTORY_STATUS_READY = 1
uint8 TRAJECTORY_STATUS_COMPLETED = 3
uint8 TRAJECTROY_STATUS_ABORT = 4
uint8 TRAJECTORY_STATUS_ILLEGAL_START = 5
uint8 TRAJECTORY_STATUS_ILLEGAL_FINAL = 6
uint8 TRAJECTORY_STATUS_IMPOSSIBLE = 7
uint8 ACTION_STOP = 8

# Its ID number will start from 1, allowing you comparing it with 0.
uint8 trajectory_flag

#zyh added
float64 pos_norm
float64 vel_norm
float64 acc_norm

# geometry_msgs/Vector3 pos_rmse
# geometry_msgs/Vector3 vel_rmse
# geometry_msgs/Vector3 acc_rmse

float64 pos_rmse_xy
float64 pos_rmse_all
float64 vel_rmse_all
float64 acc_rmse_all
