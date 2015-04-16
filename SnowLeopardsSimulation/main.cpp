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
    sc.runsimulation(//{"/Users/student/Documents/ForSim10Apr--SexF--AverageSpeed.txt","/Users/student/Documents/ForSim10Apr--SexM--AverageSpeed.txt"},
                     {"/Users/student/Documents/ForSim10Apr--SexF--MarkovTransitionsTRUE.txt","/Users/student/Documents/ForSim10Apr--SexM--MarkovTransitionsTRUE.txt"},
                     (100*100), 100, // iterations
                     "/Users/student/Documents/AverageSpeed",
                     500, {M_PI}, {10},
                     432,//432
                     0, //save move
                     15000,
                     0.00000002,//density
                     0);
    
    
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
