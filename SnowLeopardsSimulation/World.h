//
//  World.h
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 08/04/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#ifndef __SnowLeopardsSimulation__World__
#define __SnowLeopardsSimulation__World__

#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

//Header files
#include "Animal.h"
#include "Sensor.h"

class World{
    
private:

public:
    std::vector<Sensor*> trappinggrid(int, double, std::ofstream&, std::vector<double>, std::vector<double>);
    void AnimalMovement(int id, double randomstart,
                        int Iteration, int SaveMovement,
                        double CentreHome_r,double MaximumDistance,
                        double no_of_move_states,std::vector<double> probability,
                        std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
                        int NoSteps, std::ofstream &Movement ,  std::vector<Sensor*> AllSensors , std::ofstream &Captures,
                        double, int);
    double StartLocation(int, double cam_interval, double seed, double CentreHome_r);
    std::vector<std::ofstream> createfiles(std::string);
    std::string make_filename( const std::string& directory ,const std::string& basename);
    
    void oneiteration (double cam_interval,
                      int NoSteps,
                      int Iteration, int SaveMovement,
                      double CentreHome_r,double MaximumDistance,
                       double no_of_move_states,std::vector<double> probability,
                      std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
                       double seed,
                       std::ofstream &Captures, std::ofstream &Movement, std::vector<Sensor*> AllSensors, int NoAnimals  );
    void MultipleIterations(int NoSensors, int NoInterations, std::string savevalue,
                                   double cam_interval, std::vector<double> SensorWidth, std::vector<double> SensorRadius,
                                   int NoSteps,
                                   int SaveMovement,
                                   double CentreHome_r,double MaximumDistance,
                                   double no_of_move_states,std::vector<double> probability,
                                   std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector
                            ,   int NoAnimals 
                            );
};

#endif /* defined(__SnowLeopardsSimulation__World__) */
