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

// declare a random engine to be used across multiple and various method calls
static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	// x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // cout << "ParticleFilter::init" << endl;

    // This line creates a normal (Gaussian) distribution for x, y and theta
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);

    for(int a = 0; a < num_particles; a++){
        particles.push_back(Particle{
                a,
                dist_x(gen),
                dist_y(gen),
                dist_theta(gen),
                1.0

        });
    }

    // initialization is over
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find normal_distribution and default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    cout << "ParticleFilter::prediction" << endl;
    // cout << "delta_t: " << delta_t << " std_pos: [ " << std_pos[0] << " , " << std_pos[1] << " , "<< std_pos[2] << " ] velocity : " << velocity << " yawrate: " << yaw_rate << endl;


    // This line creates a normal (Gaussian) distribution for x, y and theta
    // used to add noise
    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0, std_pos[2]);


    for(int a = 0; a < num_particles; a++){
        // going straight or roughly straight
        // cout << "PRIOR Particle: [ id: " << particles[a].id << ", x: " << particles[a].x << ", y: " << particles[a].y << ", weight: " << particles[a].weight << " ]" << endl;
        if(fabs(yaw_rate) < threshold){
            double vxdt = velocity * delta_t;
            particles[a].x += vxdt * cos(particles[a].theta);
            particles[a].y += vxdt * sin(particles[a].theta);
        }
            // definitely turning
        else{
            double vdivyr = velocity / yaw_rate;
            double yrxdt = yaw_rate * delta_t; // change in theta
            particles[a].x += vdivyr * (sin(particles[a].theta + yrxdt) - sin(particles[a].theta));
            particles[a].y += vdivyr * (cos(particles[a].theta) - cos(particles[a].theta + yrxdt));
            particles[a].theta += yrxdt;
        }

        // add noise
        particles[a].x += dist_x(gen);
        particles[a].y += dist_y(gen);
        particles[a].theta += dist_theta(gen);

        // cout << "UPDATED Particle: [ id: " << particles[a].id << ", x: " << particles[a].x << ", y: " << particles[a].y << ", weight: " << particles[a].weight << " ]\n" << endl;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> landmarks, std::vector<LandmarkObs> &observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the
	// observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    // cout << "ParticleFilter::dataAssociation" << endl;

    // go through each landmark and find the one observation it is closest to
    for(int a = 0; a < landmarks.size(); a++){
        int index = 0;
        double min_distance = distance(observations[0], landmarks[a]);
        for(int b = 1; b < observations.size(); b++){
            double curr_distance = distance(observations[b], landmarks[a]);
            if(curr_distance < min_distance){
                // cout << "\t observations[b].id" << observations[b].id  << endl;

                // if this observation was already claimed for a landmark, just skip it
                if(observations[b].id > -1) continue;
                min_distance = curr_distance;
                // using this to save the index of the observation this landmark is closest to
                index = b;
            }
        }

        // save the observations closest landmark
        observations[index].id = a;
    }

    /*
    for(int a = 0; a < observations.size(); a++){
        cout << "observations[a].id: "  << observations[a].id << endl;
    }
     */
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	// more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    cout << "ParticleFilter::updateWeights" << endl;

    // used to hold transformed observations and inRange landmarks
    vector<LandmarkObs> inRange;
    vector<LandmarkObs> transformedObservations;

    // shared values in weights calculation
    const double sigxp2x2 = 2 * pow(std_landmark[0], 2);
    const double sigyp2x2 = 2 * pow(std_landmark[1], 2);
    const double pix2xsigxxsigy = 2 * M_PI * std_landmark[0] * std_landmark[1];

    // for each particle
    for(int a = 0; a < num_particles; a++){
        // cout << "a: " << a << endl;
        // find the landmarks in sensor range
        findLandmarksinRange(particles[a], map_landmarks, sensor_range, inRange);

        // transform observations to map coordinate space
        transformToCoordinateSpace(particles[a], observations, transformedObservations);

        // associate observations with landmarks
        dataAssociation(inRange, transformedObservations);

        // update weights
        for(int b = 0; b < transformedObservations.size(); b++){
            // cout << "id: " << transformedObservations[b].id << endl;
            // skip unaffliated
            if(transformedObservations[b].id == -1) continue;

            double diff_x = transformedObservations[b].x - inRange[transformedObservations[b].id].x;
            // cout << "diff_x: " << diff_x << endl;
            double diff_y = transformedObservations[b].y - inRange[transformedObservations[b].id].y;
            // cout << "diff_y: " << diff_y << endl;
            double power = -1 * (pow(diff_x, 2) / sigxp2x2 + pow(diff_y, 2) / sigyp2x2);
            // cout << "power: " << power << endl;
            double e = exp(power);
            // cout << "e: " << e << endl;
            double weight = e / pix2xsigxxsigy;
            // cout << "weight: " << weight << endl;
            if (weight > 0) particles[a].weight *= weight; // TODO see if this is necessary
        }
        cout <<"\tparticle[" << a << "]'s weight: " << particles[a].weight << endl;
        if(particles[a].weight > max_weight) max_weight = particles[a].weight;
    }

    // reset the max weight to prevent infinite loops in the future
    max_weight = 0;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    cout << "ParticleFilter::resample" << endl;

    // TODO find source of infinite loop
    // return;

    vector<Particle> resampled;

    // generate random numbers
    std::default_random_engine gen;

    std::uniform_real_distribution<double> distribution(0.0, 2 * max_weight);
    std::uniform_int_distribution<> range(0, num_particles-1);

    int index = range(gen);
    double beta = distribution(gen);

    // resample until i hae enough particles
    for (int i = 0; i < num_particles; i++) {
        beta += distribution(gen);
        while (beta > particles[index].weight) {
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }
        cout << "\tparticle[" << index << "]'s weight " << particles[index].weight << endl;
        resampled.push_back(particles[index]);
    }

    particles = resampled;
}

double ParticleFilter::distance(LandmarkObs observed, LandmarkObs actual){
    return sqrt(pow(observed.x - actual.x, 2) + pow(observed.y - actual.y, 2));
}

void ParticleFilter::findLandmarksinRange(Particle particle, Map landmarks, double sensor_range, std::vector<LandmarkObs> &inRange){
    for(int a = 0; a < landmarks.landmark_list.size(); a++){
        double distance = sqrt(pow(particle.x - landmarks.landmark_list[a].x_f, 2) + pow(particle.y - landmarks.landmark_list[a].y_f, 2));

        if(distance <= sensor_range){
            inRange.push_back(LandmarkObs{
                    landmarks.landmark_list[a].id_i,
                    landmarks.landmark_list[a].x_f,
                    landmarks.landmark_list[a].y_f
            });
        }
    }
}

void ParticleFilter::transformToCoordinateSpace(Particle particle, std::vector<LandmarkObs> observations, std::vector<LandmarkObs> &transformedObservations){
    // cout << "ParticleFilter::transformToCoordinateSpace" << endl;
    for(int a = 0; a < observations.size(); a++){
        transformedObservations.push_back(LandmarkObs{
                -1,
                particle.x + cos(particle.theta) * observations[a].x - sin(particle.theta) * observations[a].y,
                particle.y + sin(particle.theta) * observations[a].x + cos(particle.theta) * observations[a].y
        });
        // cout << transformedObservations[a].id << endl;
    }
}

// do not mess with the following
Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y) {
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

std::string ParticleFilter::getAssociations(Particle best) {
    std::vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

std::string ParticleFilter::getSenseX(Particle best) {
    std::vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    std::string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

std::string ParticleFilter::getSenseY(Particle best) {
    std::vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    std::string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
