/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// Assign weights
	num_particles = 200;
	double init_weight = 1;
	weights.assign(num_particles, init_weight);

	// Initialize all particles to first position, with some noise
	std::default_random_engine gen;
	double std_x, std_y, std_psi;
	std_x = std[0];
	std_y = std[1];
	std_psi = std[2];

	std::normal_distribution<double> dist_x(x,std_x);
	std::normal_distribution<double> dist_y(y,std_y);
	std::normal_distribution<double> dist_psi(theta,std_psi);

	Particle temp_p;

	for(int i = 0; i < num_particles; i++){
		temp_p.id = i;
		temp_p.x = dist_x(gen);
		temp_p.y = dist_y(gen);
		temp_p.theta = dist_psi(gen);
		temp_p.weight = 1.0;
		// Push back
		particles.push_back(temp_p);
	}

	is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;
	for(int i = 0; i < particles.size(); i++){
		double x_temp = particles[i].x;
		double y_temp = particles[i].y;
		double theta_temp = particles[i].theta;
		// Split into two cases
		if(fabs(yaw_rate) > 0.001){
			x_temp = x_temp + velocity/yaw_rate * ( sin (theta_temp + yaw_rate*delta_t) - sin(theta_temp));
        	y_temp = y_temp + velocity/yaw_rate * ( cos(theta_temp) - cos(theta_temp+yaw_rate*delta_t) );
        	theta_temp = theta_temp + yaw_rate*delta_t;
		}
		else{
			x_temp = x_temp + velocity*delta_t*cos(theta_temp);
			y_temp = y_temp + velocity*delta_t*cos(theta_temp);
		}
		// Add some noise to measurement
		std::normal_distribution<double> dist_x(x_temp,std_pos[0]);
		std::normal_distribution<double> dist_y(y_temp,std_pos[1]);
		std::normal_distribution<double> dist_psi(theta_temp,std_pos[2]);
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_psi(gen);
		//std::cout << "p.x = " << particles[i].x << " | " << "p.theta = " << particles[i].theta << " | ";
		//std::cout << std::endl;
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Not implemented
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	// clear previous computed weights
	weights.clear();
	// define standard deviations
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];
	double std_x_sq2 = std_x*std_x;
	double std_y_sq2 = std_y*std_y;
	long double constant_1 = 1 / (2 * M_PI * std_x * std_y);

	// for each particle
	for(int i = 0; i < particles.size(); i++){
		// save particle information for later use
		double theta_temp = particles[i].theta;
		double x_temp = particles[i].x;
		double y_temp = particles[i].y;
		long double weight_p = 1.0;
		// given a particle's location, find all transformed measurement
		for(int j = 0; j < observations.size(); j++){
			// measurement in vehicle coordinates
			double x_v = observations[j].x;
			double y_v = observations[j].y;

			// coordinate transformation starts here
			// 1. Rotation
			double x_m = x_v*cos(theta_temp) - y_v*sin(theta_temp);
			double y_m = x_v*sin(theta_temp) + y_v*cos(theta_temp);
			// 2. Add translational vector
			x_m += x_temp;
			y_m += y_temp;
			// See if the transformed measuremnt is within the sensor range.
			// if not within sensor range, skip
			if(dist(x_temp, y_temp, x_m, y_m) > sensor_range)
				continue;

			// Setup a tempory minimum id for this observation:
			int min_id = -1;
			double threshhold = 50;
			// given a transformed measurement find the closest landmark
			for(int k = 0; k < map_landmarks.landmark_list.size(); k++){
				// save each landmark's x,y coordiante for computing distance
				double x_land = map_landmarks.landmark_list[k].x_f;
				double y_land = map_landmarks.landmark_list[k].y_f;
				int id_land = map_landmarks.landmark_list[k].id_i-1;
				double smallest_distance = dist(x_m, y_m, x_land, y_land);
				if (smallest_distance < threshhold){
					min_id = id_land;

					threshhold = smallest_distance;
				}
			}
			//std::cout << "min id : " << min_id << std::endl;
			if(min_id >= 0){
				// Update the weight
				long double u_x = map_landmarks.landmark_list[min_id].x_f;
				long double u_y = map_landmarks.landmark_list[min_id].y_f;
				long double x_diff_sq2 = (x_m - u_x)*(x_m - u_x);
				long double y_diff_sq2 = (y_m - u_y)*(y_m - u_y);
				long double part2 = x_diff_sq2 / std_x_sq2 + y_diff_sq2 / std_y_sq2;
				long double weight_temp = constant_1*exp( (-1/2.) * (part2));
				//std::cout << "weight temp : " << weight_temp << std::endl;
				weight_p *= weight_temp;
			}
		}
		particles[i].weight = weight_p;
		//std::cout << "#" << i << ": " << weight_p << std::endl;
		weights.push_back(weight_p);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine gen;
	std::discrete_distribution<int> distribution(weights.begin(), weights.end());
	int weighted_index = distribution(gen);

	std::vector<Particle> resampled_particles;

	for(int i = 0; i < particles.size(); i++){
		resampled_particles.push_back(particles[weighted_index]);
	}
	particles = resampled_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
