//
//  main.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 25/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include <iostream>
#include <random>
#include "Animal.h"
#include "Sensor.h"
#include "World.h"


int main(int argc, const char * argv[]) {
    // insert code here...

    World sc;
    //"TestReadIn.txt"
    sc.runsimulation({"/Users/student/Documents/ForSim30Apr--SexF--MarkovTransitionsTRUE.txt","/Users/student/Documents/ForSim30Apr--SexM--MarkovTransitionsTRUE.txt"},
                     (1*1), 1, "/Users/student/Documents/xx",
                     500, {M_PI}, {200},
                     432, //432
                     1, //save move
                     15000,
                     0.000002,//10000,
                     0.5);
    
    
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
