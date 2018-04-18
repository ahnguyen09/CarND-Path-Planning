#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

#define SCAN_TIME 0.02 // define scan time of 20ms
#define MAX_SPEED 22.12848 // m/s -> limit max speed to 49.5 mph
#define MAX_ACCEL 3.0 // set max accel to 5.0 m/s/s
#define JERK 3.0 // set jerk rate to 5.0 m/s/s/s
#define MPH_TO_MS 0.44704 // multipler to convert MPH to MS

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

	//start in lane 1
	int lane = 1;

	//starting reference velocity
	double ref_vel = 0.0; //m/s -> start at 0 mph

	//starting reference accel
	double ref_accel = 0.0; //m/s/s

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&ref_vel,&ref_accel]
							(uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"]; // in terms of mph

						//cout << car_speed << endl;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
						// sensor_fusion[i][5] - [ id, x, y, vx, vy, s, d]

          	json msgJson;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
						int prev_size = previous_path_x.size();

						// if previous path exist, change car s to previous ending waypoint s
            if (prev_size > 0) {
              car_s = end_path_s;
            }

						
						//create widely spaced point so we can interpolate between the point to move the car smoothly
						vector<double> ptsx;
            vector<double> ptsy;

						//reference will be either current current car's orientation or previous path orientation
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            // if previous size is almost empty usre car orientaiton as reference
            if ( prev_size < 2 ) {
                //assume previous point is along the saw yaw
                double prev_car_x = ref_x - cos(ref_yaw);
                double prev_car_y = ref_y - sin(ref_yaw);

                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);

                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);
            } else {
                // reference orientation is a deepest waypoint in the que
								// this overrites reference state
                ref_x = previous_path_x[prev_size - 1];
                ref_y = previous_path_y[prev_size - 1];

                double ref_x_prev = previous_path_x[prev_size - 2];
                double ref_y_prev = previous_path_y[prev_size - 2];
                ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);

                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);
            }

            // Setting up target points in ahead in Frenet
            vector<double> next_wp0 = getXY(car_s + 15, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s + 30, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s + 45, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            // transforming from global x,y to local car coordinates at each waypoint
            for ( int i = 0; i < ptsx.size(); i++ ) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = shift_x*cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw);
              ptsy[i] = shift_x*sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw);
            }

            // Create the spline.
            tk::spline s;
            s.set_points(ptsx, ptsy);

            // points to send to simulator
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

						// start from the points that was already calculated
            for ( int i = 0; i < prev_size; i++ ) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

						// Set Target position to be 30m in front of us when there are no cars ahead
						double target_x = 30.0;
						double target_speed = 22.12848; //m/s - > set target speed to 49.5 mph
						double const_time_gap = 15;

						float closest_s_diff = 99999;
						float closest_s_id = 0;

						//sensor fusion to determine target position and speed
            for ( int i = 0; i < sensor_fusion.size(); i++ ) {
							float s = sensor_fusion[i][5];
              float d = sensor_fusion[i][6];

							if ((d > 4*lane && d < 4*(lane+1)) && abs(s-car_s) <= 300) { // same lane and within 300m
								// Find car speed.
								double vx = sensor_fusion[i][3];
								double vy = sensor_fusion[i][4];
								double check_speed = sqrt(vx*vx + vy*vy);
								double check_car_s = sensor_fusion[i][5];
								// Estimate car s position after executing previous trajectory.
								check_car_s += ((double)prev_size*SCAN_TIME*check_speed);

							
								car_speed = car_speed*MPH_TO_MS;
								// change target x and target speed based on constant time gap and sensor fusion 

								// save closes id of car that is ahead
								if ((check_car_s-car_s) < closest_s_diff && check_car_s > car_s) {
									closest_s_diff = check_car_s-car_s;
									closest_s_id = i;
								}

								if (check_car_s > car_s && (check_car_s-car_s) < 30 && i == closest_s_id) {
									const_time_gap = 5 + car_speed*(3.5/22.352) - (check_speed-car_speed)*(0.5/22.352);
									target_x = check_car_s - car_s - const_time_gap;
									
									target_speed = check_speed - 0.2;
								}
							}               
            }

            double target_y = s(target_x); //evaluate spine at target x
            double target_dist = sqrt(target_x*target_x + target_y*target_y);

            double x_add_on = 0;

						//increment x in short increments based on reference velocity
            for( int i = 1; i < 50 - prev_size; i++ ) {
							//vf = vi + a*t + 0.5*j*t^2
							//af = ai + j*t

							ref_vel += ref_accel*SCAN_TIME + 0.5*JERK*SCAN_TIME*SCAN_TIME;
							//limit velocity
							if ( ref_vel > MAX_SPEED ) {
                ref_vel = MAX_SPEED;
							}

							/*
								target speed > ref_vel, target speed < ref_vel & target speed - ref_vel < 0.5*ref_accel^2/jerk
									+jerk
								target speed < ref_vel, target speed > ref_vel & target speed - ref_vel < 0.5*ref_accel^2/jerk
									-jerk
								== 
									accel = 0;
							*/ 

							if (target_speed > ref_vel || 
							    (target_speed < ref_vel && abs(target_speed-ref_vel) <= abs(0.5*ref_accel*ref_accel/JERK)) ) {
								ref_accel += JERK*SCAN_TIME;
							}
							else if (target_speed < ref_vel || 
											 (target_speed > ref_vel && abs(target_speed-ref_vel) <= abs(0.5*ref_accel*ref_accel/JERK)) || 
											 (closest_s_diff < const_time_gap)  ) {
								ref_accel -= JERK*SCAN_TIME;
							}
							/*
							else {
								ref_accel = 0;
							}
							*/

							//limit accel
							if ( ref_accel > MAX_ACCEL ) {
                ref_accel = MAX_ACCEL;
							}

							/*
              ref_vel += speed_diff;
              if ( ref_vel > MAX_SPEED ) {
                ref_vel = MAX_SPEED;
              } else if ( ref_vel < MAX_ACC ) {
                ref_vel = MAX_ACC;
              }
							*/

							/*
							cout << "target x " << target_x << endl;
							
							cout << "target speed " << target_speed << endl; 

							cout << "ref speed " << ref_vel << endl;
							cout << "ref accel " << ref_accel << endl;
							*/

              double N = target_dist/(SCAN_TIME*ref_vel); //distance travel / distance travel at reference speed
              double x_point = x_add_on + target_x/N;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

							//rotate local car coordinate back to world coordinate
              x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
              y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

						msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
