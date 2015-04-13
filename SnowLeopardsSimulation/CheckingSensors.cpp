//
//  CheckingSensors.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 09/04/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include "CheckingSensors.h"
#include "Sensor.h"

#include <random>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>

CheckingSensor::CheckingSensor(){};
CheckingSensor::CheckingSensor(double currentx,double currenty, double currentangle, double previousx,double previousy,
                         std::vector<Sensor*> AllSensors,
                         std::ofstream &Movement,
                         std::ofstream &Captures,
                         int stepcounter, int id, int iterationnumber,
                         int NoSensors){
    
    
    for(int sensor=0; sensor<(NoSensors); sensor++){
        double x = AllSensors[sensor] -> get_x();double y = AllSensors[sensor] -> get_y();
        if(sqrt(pow(x-previousx,2)+pow(y-previousy,2))<15000){
        AllSensors[sensor] -> CapturesIntersection(currentx,currenty,
                                                   previousx,previousy,
                                                   id,currentangle,
                                                   iterationnumber, Captures,
                                                   stepcounter);
        }
    }; // END of sensor for loop
};