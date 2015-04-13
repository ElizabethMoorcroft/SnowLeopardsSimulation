//
//  CheckingSensors.h
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 09/04/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#ifndef __SnowLeopardsSimulation__CheckingSensors__
#define __SnowLeopardsSimulation__CheckingSensors__

#include "Sensor.h"
#include <stdio.h>
#include <random>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>

class CheckingSensor{
public:
    CheckingSensor();
    CheckingSensor(double currentx,double currenty, double currentangle, double previousx,double previousy,
                   std::vector<Sensor*> AllSensors,
                             std::ofstream &Movement,
                             std::ofstream &Captures,
                             int stepcounter, int i, int iterationnumber,
                   int);
};

#endif /* defined(__SnowLeopardsSimulation__CheckingSensors__) */
