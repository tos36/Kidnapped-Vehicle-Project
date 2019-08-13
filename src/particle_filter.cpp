/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // std::cout << "init start" <<std::endl;
  
  num_particles = 100;  // TODO: Set the number of particles
  weights.resize(num_particles);
  
  // Gaussian distribution for x, y, theta
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  // vector<Particle> list_particles; //list of particles
  for(int i=0; i < num_particles; ++i){
  	Particle temp;
  	temp.id = i;
    temp.x = dist_x(gen);
    temp.y = dist_y(gen);
    temp.theta = dist_theta(gen);
    temp.weight = 1.0;
    particles.push_back(temp);
    weights[i] = 1.0;
    
  }
  is_initialized = true;

  // std::cout << "init end" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // std::cout << "pred start" << std::endl;
  
  
  // Gaussian distribution for x, y, theta
  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  for(int i=0; i < num_particles; ++i){
    double theta = particles[i].theta;
    
    if(fabs(yaw_rate)<0.0001){
      // avoid zero-division
      particles[i].x += velocity * delta_t * cos(theta);
      particles[i].y += velocity * delta_t * sin(theta);
    }else{
      particles[i].x += velocity * (sin(theta + yaw_rate * delta_t) - sin(theta)) / yaw_rate;
      particles[i].y += velocity * (cos(theta) - cos(theta + yaw_rate * delta_t)) / yaw_rate;
      particles[i].theta += yaw_rate * delta_t;
    }
    
    // Adding noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    
  }

  // std::cout << "pred end" << std::endl;
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // std::cout << "DA start" << std::endl;
  double dist_temp;
  unsigned int num_obs = observations.size();
  unsigned int num_pred = predicted.size();
  
  for(int i=0; i<num_obs; ++i){
    double dist_min = 99999;
    int nearest_id = -1;
    for(int j=0; j<num_pred; ++j){
      dist_temp = dist(observations[i].x, observations[i].y,
                       predicted[j].x, predicted[j].y);
      
      if (dist_temp < dist_min){
        dist_min = dist_temp;
        nearest_id = predicted[j].id;
      }
    }
    observations[i].id = nearest_id;
  }

  // std::cout << "DA end" << std::endl;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // std::cout << "update start" << std::endl;
  for (int i=0; i<num_particles; ++i){
    vector<LandmarkObs> transformed_observations;
    Particle& p = particles[i];
    p.weight = 1.0;
    
    // transform coordinate
    for (int j=0; j<observations.size(); ++j){
      LandmarkObs obs = observations[j];
      LandmarkObs obs_trans;
      obs_trans.x = obs.x * cos(p.theta) - obs.y * sin(p.theta) + p.x;
      obs_trans.y = obs.x * sin(p.theta) + obs.y * cos(p.theta) + p.y;
      obs_trans.id = obs.id;
      transformed_observations.push_back(obs_trans);
    }
    
    // pick up landmarks in senser range
    vector<LandmarkObs> predicted_landmarks;
    
    for (int j=0; j<map_landmarks.landmark_list.size(); j++){
      double distance;
      double land_x = map_landmarks.landmark_list[j].x_f;
      double land_y = map_landmarks.landmark_list[j].y_f;
      distance = dist(p.x, p.y, land_x, land_y);
      
      if (distance < sensor_range){
        LandmarkObs land_pred;
        land_pred.id = map_landmarks.landmark_list[j].id_i;
        land_pred.x = map_landmarks.landmark_list[j].x_f;
        land_pred.y = map_landmarks.landmark_list[j].y_f;
        predicted_landmarks.push_back(land_pred);
      }
    }
    
    // Associate observations to prediction
    dataAssociation(predicted_landmarks, transformed_observations);
    
    double prob = 1.0;
    for (int j=0; j<predicted_landmarks.size(); ++j){
      double dist_min = sensor_range;
      int id_min = -1;
      
      double pred_x = predicted_landmarks[j].x;
      double pred_y = predicted_landmarks[j].y;
      
      for (int k=0; k<transformed_observations.size(); k++){
        double trans_x = transformed_observations[k].x;
        double trans_y = transformed_observations[k].y;
        double dist_temp = dist(pred_x, pred_y, trans_x, trans_y);
        if (dist_temp < dist_min){
          dist_min = dist_temp;
          id_min = k;
        }
      }
      
      double std_x = std_landmark[0];
      double std_y = std_landmark[1];
      double c = 1/(2 * M_PI * std_x * std_y);
      
      if (id_min != -1){
        double trans_x = transformed_observations[id_min].x;
        double trans_y = transformed_observations[id_min].y;
        double multi = c * exp(-0.5 * (pow((trans_x - pred_x)/std_x, 2) 
                                       + pow((trans_y - pred_y)/std_y, 2)));   
        prob *= multi;
    }   
  }
     p.weight = prob;
     weights[i] = prob;
  }

  // std::cout << "update end" << std::endl;
}




void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   // std::cout << weights << std::endl;
  
  
  vector<Particle> new_particles;
  default_random_engine gen;
  std::discrete_distribution<int> index(weights.begin(), weights.end());
  
  for (int i=0; i<num_particles; ++i){
    new_particles.push_back(particles[index(gen)]);
  }
  particles = new_particles;
   

 //  std::cout << "resample end" << std::endl;
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