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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	//Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	default_random_engine gen;

	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	std_x = std[0];	std_y = std[1]; std_theta = std[2]; //set the standard deviations

	// This line creates a normal(Gaussian) distributions for x, y and z.
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);


	for (unsigned int i = 0; i < num_particles; ++i) {
	//Sample from the x,y and z normal distrubtions and concatenate to create particles

		Particle P = {};
		P.id = i;
		P.x = dist_x(gen);
		P.y = dist_y(gen);
		P.theta = dist_theta(gen);

		//push particle to vector
		particles.push_back(P); 
	}

	//set initialization flag
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.

	default_random_engine gen;

	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	std_x = std_pos[0];	std_y = std_pos[1]; std_theta = std_pos[2]; //set the standard deviations
	
	double v = velocity; 
	double thetad = yaw_rate; 

	
	for (unsigned int i = 0; i < num_particles; ++i) {
		//itterate through and predict particles based on motion model

		//current position and orientation of particle
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		//updates position and orientation of particle
		double new_x, new_y, new_theta;

		if (thetad == 0) {
			//for zero rotation rate
			new_x = x + v * delta_t * cos(theta);
			new_y = y + v * delta_t * sin(theta);
			new_theta = theta;
		}
		else {
			//for nonzero rotation rate
			new_x = x + (v / thetad) * (sin(theta + thetad * delta_t) - sin(theta));
			new_y = y - (v / thetad) * (cos(theta + thetad * delta_t) - cos(theta));;
			new_theta = theta + thetad * delta_t;
		}

		//create norma;l distributions of the new positions and orientation
		normal_distribution<double> dist_x(new_x, std_x);
		normal_distribution<double> dist_y(new_y, std_y);
		normal_distribution<double> dist_theta(new_theta, std_theta);

		//update particles based on samples from the normal distributions
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		
	}

}

void ParticleFilter::dataAssociation(const Map map_landmarks, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: both observations and landmarks should be in the same coordinate frame



	for (unsigned int i = 0; i < observations.size(); i++) {
		//iterate through observations

		//initialize distance to closest landmark as a large value
		double shortest_distance = 1000000000;

		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			//iterate through landmarks

			//calculate distance between landmark and observation
			double distance = sqrt(pow(map_landmarks.landmark_list[j].x_f - observations[i].x, 2) + pow(map_landmarks.landmark_list[j].y_f - observations[i].y, 2));

			if (distance < shortest_distance) {
				//if new closer landmark is found

				//update the shortest distance variable
				shortest_distance = distance;
				//set the id of the observation as the order of the closest landmark
				observations[i].id = j;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//clear weights vector
	ParticleFilter::weights.clear();


	for (unsigned int i = 0; i < num_particles; ++i) {
		//iterate through particles

		//particle coordinates
		double xp = particles[i].x;
		double yp = particles[i].y;
		double thetap = particles[i].theta;

		//Vector of observations transformed to global coordinates
		vector<LandmarkObs> map_cor_observations;

		for (unsigned int j = 0; j < observations.size(); ++j) {
			//iterate through observations

			//define landmark object
			LandmarkObs map_cor_observation; 

			//Apply coordinate transformation from vehicle coordinates to global coordinates
			map_cor_observation.id = observations[j].id;
			map_cor_observation.x = xp + cos(thetap)*observations[j].x - sin(thetap)*observations[j].y;
			map_cor_observation.y = yp + sin(thetap)*observations[j].x + cos(thetap)*observations[j].y;

			//push transformed object to vector
			map_cor_observations.push_back(map_cor_observation);
		}

		//Associate each observation to the closest landmark
		dataAssociation(map_landmarks, map_cor_observations);

		//initialize particle weight as unity
		double particle_weight = 1;

		//standard deviations of landmark positions
		double std_land_x = std_landmark[0];
		double std_land_y = std_landmark[1];

		//vectors for the association of each particle
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		for (unsigned int j = 0; j < observations.size(); ++j) {
			//iterate through observations

			//global coordinates of the landmark associated with each observation
			double predicted_x = map_landmarks.landmark_list[map_cor_observations[j].id].x_f;
			double predicted_y = map_landmarks.landmark_list[map_cor_observations[j].id].y_f;

			//global coordinates of observation
			double observe_x = map_cor_observations[j].x;
			double observe_y = map_cor_observations[j].y;
			
			//calculate particles weight
			particle_weight *= (1 / (2 * M_PI*std_land_x*std_land_y)) * exp(-(pow(predicted_x - observe_x, 2) / (2 * pow(std_land_x, 2))) - (pow(predicted_y - observe_y, 2) / (2 * pow(std_land_y, 2))));
		

			//update associations vector
			associations.push_back(map_cor_observations[j].id + 1);
			sense_x.push_back(observe_x);
			sense_y.push_back(observe_y);
		}

		//update particle weight
		particles[i].weight = particle_weight;
		weights.push_back(particle_weight);

		//relate the associations vector to particle
		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);

	}


}

void ParticleFilter::resample() {
	//Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//Define vector for resampled particles
	vector<Particle> resampled_particles;

	//random generator & distribution
	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> disc_dist(weights.begin(), weights.end());

	for (unsigned int i = 0; i < num_particles; i++) {
		//obtain and push random samples 
		resampled_particles.push_back(particles[disc_dist(gen)]);
	}

	//update particles vector
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, const std::vector<int> associations, 
                                     const std::vector<double> sense_x, const std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	Particle new_particle = particle;
	new_particle.associations= associations;
	new_particle.sense_x = sense_x;
	new_particle.sense_y = sense_y;

	return new_particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
