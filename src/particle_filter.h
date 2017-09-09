/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

using namespace std;

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	vector<int> associations;
	vector<double> sense_x;
	vector<double> sense_y;
};



class ParticleFilter {
private:
	// Number of particles to draw
	int num_particles;
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	vector<double> weights;

	// threshold value
	float threshold;

	/**
	 * Method to get distance between 2D points
	 *
	 * @param observed
	 * @param actual
	 * @return distance between two points
	 */
	double distance(LandmarkObs observed, LandmarkObs actual);

	/**
	 * Method to find landmarks in sensor range
	 *
	 * @param particle
	 * @param landmarks
	 * @param sensor_range
	 * @param inRange
	 */
	void findLandmarksinRange(Particle particle, Map landmarks, double sensor_range, vector<LandmarkObs> &inRange);

	/**
	 * Method to transform observations-in-relation-to-car to map coordinates
	 *
	 * @param particle
	 * @param observations
	 * @param transformedObservations
	 */
	void transformToCoordinateSpace(Particle particle, vector<LandmarkObs> observations, vector<LandmarkObs> &transformedObservations);
	
public:
	
	// Set of current particles
	vector<Particle> particles;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter() : num_particles(100), is_initialized(false), threshold(0.001) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param landmarks Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(vector<LandmarkObs> landmarks, vector<LandmarkObs> &observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], vector<LandmarkObs> observations,
			Map map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, vector<int> associations, vector<double> sense_x, vector<double> sense_y);
	
	string getAssociations(Particle best);
	string getSenseX(Particle best);
	string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}
};



#endif /* PARTICLE_FILTER_H_ */
