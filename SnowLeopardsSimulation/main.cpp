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
<<<<<<< HEAD
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

=======
    World sc;
    //"TestReadIn.txt"
    sc.runsimulation({"/Users/student/Documents/ForSim10Apr--SexF--MarkovTransitionsTRUE.txt","/Users/student/Documents/ForSim10Apr--SexM--MarkovTransitionsTRUE.txt"},
                     (1*1), 100, "/Users/student/Documents/xx",
                     500, {M_PI}, {200},
                     432, //432
                     0, //save move
                     15000,
                     100000,//10000,
                     0.5);
    
    
>>>>>>> parent of 9f11446... Working correct code
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
