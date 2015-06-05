//
//  Animal.h
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 25/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#ifndef __SnowLeopardsSimulation__Animal__
#define __SnowLeopardsSimulation__Animal__

#include <stdio.h>
#include <vector>
#include "Sensor.h"

class Animal{
private:
    // Identifiers
    int AnimalId;
    int Iteration;
    int step_number;
    int SaveMovement;
    int Sex;
    
    // movement setups
    double CentreHome_x;
    double CentreHome_y;
    double CentreHome_r;
    double MaximumDistance;
    double no_of_move_states;
    std::vector<double> probability;
    std::vector<std::vector<double>> mean_vector;
    std::vector<std::vector<double>> variance_vector;
    std::vector<std::vector<double>> transitions_vector;
    
    // current movement
    int Current_state;
    double Current_x;
    double Current_y;
    double Current_angle;
    double Current_distance;

    
    //overall movement
    double Total_distance;

public:
    
    Animal();
    Animal(int AnimalId,int Iteration, int SaveMovement,
                   // movement setups
           double CentreHome_r,double MaximumDistance,double no_of_move_states,
           std::vector<double> probability, std::vector<std::vector<double>> mean_vector, std::vector<std::vector<double>> variance_vector,
           std::vector<std::vector<double>> transitions_vector,
           double x, double y, double a, double seed,
           std::ofstream &Movement, std::vector<Sensor*> AllSensors , std::ofstream &Captures, int);
    
    void set_AnimalId(int Id){AnimalId= Id;};
    void set_MaximumDistance(double newmax){MaximumDistance = newmax;};
    
    void set_Current_x(double newx){Current_x = newx;};
    void set_Current_y(double newy){Current_y = newy;};
    void set_Current_state(double news){Current_state = news;};
    
    void set_no_of_move_states(double nostates){no_of_move_states=nostates;};
    void set_probability(std::vector<double> probs){probability =probs;};
    void set_mean_vector(std::vector<std::vector<double>> mean){mean_vector = mean;};
    void set_variance_vector(std::vector<std::vector<double>> vect){variance_vector = vect;};
    
    void set_CentreHome_x(double newx){CentreHome_x = newx;};
    void set_CentreHome_y(double newy){CentreHome_y = newy;};
    void set_CentreHome_r(double newr){CentreHome_r = newr;};

    double get_Current_state(){return Current_state;};
    double get_Current_x(){return Current_x;};
    double get_Current_y(){return Current_y;};
    double get_Total_distance(){return Total_distance;};
    
    
    void LocationVector(double previous_x, double previous_y, int LeaveEntreCode, int End, std::ofstream &Movement, std::vector<Sensor*> AllSensors , std::ofstream &Captures);
    void UpdateLocation (double, std::ofstream &Movement, std::vector<Sensor*>, std::ofstream &Captures);
    void NewLocation (double);
    std::vector<double> NewLocationFromMult(double seed, std::vector<double> mean, std::vector<double> variance); //tested
    std::vector<double> MultivariteNorm (double, std::vector <double>, std::vector <double>); //tested
    std::vector<double> MatrixMultip (std::vector <double>, std::vector <double> ); // tested
    int PickState (double,  double, std::vector<double>); // tested
    std::vector<double> ThreeSeeds(double seed);
    double Normal (double, double, double); //test
    
    double CalNext_X(double, double); // tested
    double CalNext_Y(double, double); // tested
    double RangeAngle(double); // tested
    

};

#endif /* defined(__SnowLeopardsSimulation__Animal__) */
