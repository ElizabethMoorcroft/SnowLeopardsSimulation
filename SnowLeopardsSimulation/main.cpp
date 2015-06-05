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
    World sc; World sc1; World sc2; World sc3;
    sc.runsimulation({"/Users/student/Documents/ForSim16Apr--SexF--MarkovTransitionsTRUE.txt"},
                     (1*1), 1, // iterations
                     "/Users/student/Documents/TEST",//TransitionsTRUE
                     4000, {M_PI}, {10},
                     432,//432
                     1, //save move
                     15000,
                     0.00000002,//density
                     0.5);

    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
