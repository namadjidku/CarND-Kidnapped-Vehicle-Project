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
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  
  // This line creates a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  // Set of current particles
  for(int i = 0; i < num_particles; i++) {
    Particle particle = Particle();
    particle.x = dist_x(gen); 
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;
    particle.id = i;
  	particles.push_back(particle);
    
    // to remove
    weights.push_back(1.0);
  }
  
  is_initialized = true; 

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine   useful.
   *  http://en.cp ference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for(int i = 0; i < num_particles; i++) {
    
    if (fabs(yaw_rate) == 0.0) {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
      particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
    } 
    else {  
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + dist_x(gen); 
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta  + yaw_rate*delta_t)) + dist_y(gen); 
      particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
    }
  }
  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  
  for(int i = 0; i < observations.size(); i++) 
  {
    if(predicted.size() > 0) 
    {
      double closest_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
	  observations[i].id = 0;
      
      for(int j = 1; j < predicted.size(); j++) {

        double temp = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

        if(temp < closest_dist) {
          closest_dist = temp;
          observations[i].id = j;
        }
      }
    }
  }
}

void ParticleFilter::transformObservations(const std::vector<LandmarkObs> observations, std::vector<LandmarkObs> &map_observations, Particle particle) {

	for(int j = 0; j < observations.size(); j++) {
       // homogeneous transformation
       LandmarkObs landmarkObs = LandmarkObs();
       landmarkObs.x = particle.x  + cos(particle.theta) * observations[j].x - sin(particle.theta) * observations[j].y;
       landmarkObs.y = particle.y  + sin(particle.theta) * observations[j].x + cos(particle.theta) * observations[j].y;
       map_observations.push_back(landmarkObs);
     }
}

void ParticleFilter::findLandmarksInSensorRange(double sensor_range, const Map &map_landmarks, 
                                                  std::vector<LandmarkObs> &landmarks_in_range, Particle particle) {

    for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      
       double temp = dist(particle.x, 
                          particle.y, 
                          map_landmarks.landmark_list[j].x_f, 
                          map_landmarks.landmark_list[j].y_f);
      
      if(temp <= sensor_range) {
        LandmarkObs landmarkObs = LandmarkObs();
        landmarkObs.id =  map_landmarks.landmark_list[j].id_i;
        landmarkObs.x =  map_landmarks.landmark_list[j].x_f;
        landmarkObs.y =  map_landmarks.landmark_list[j].y_f;
        landmarks_in_range.push_back(landmarkObs);
      }
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	// Update each particle weight
	for(int i = 0; i < num_particles; i++) {

		// Transform map observations from vehicle coordinates to map coordinates
		vector<LandmarkObs> map_observations;
		transformObservations(observations, map_observations, particles[i]);

		// Find map landmarks within the sensor range
		vector<LandmarkObs> landmarks_in_range;
		findLandmarksInSensorRange(sensor_range, map_landmarks, landmarks_in_range, particles[i]);

		// Associate landmark id to each observation landmark 
      
		dataAssociation(landmarks_in_range, map_observations);
      
       //********************* FOR DEBUGING//
         
      	vector<int> associations;
     	vector<double> sense_x;
     	vector<double> sense_y;
     
     	for(int j = 0; j < map_observations.size(); j++) {
          associations.push_back(landmarks_in_range[map_observations[j].id].id);
          sense_x.push_back(map_observations[j].x);
          sense_y.push_back(map_observations[j].y);
     	}
     
     	SetAssociations(particles[i], associations, sense_x, sense_y);

        //********************* FOR DEBUGING//
      
      	// updating the weights 
         double probability = 1.0; 

         for(int j = 0; j < map_observations.size(); j++) {
           // homogeneous transformation
           probability *=  multiv_prob(std_landmark[0], std_landmark[1], 
                                       map_observations[j].x, map_observations[j].y, 
                                       landmarks_in_range[map_observations[j].id].x, landmarks_in_range[map_observations[j].id].y); 
         }

         weights[i] = particles[i].weight = probability;
	}
}

void ParticleFilter::resample() {

	/**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	std::discrete_distribution<int> resample_d(weights.begin(), weights.end());
	std::vector<Particle> new_particles;

	for(int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[resample_d(gen)]);
	}

	particles = new_particles;
}


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}