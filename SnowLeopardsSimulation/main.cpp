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
    sc.runsimulation({"/Users/student/Documents/changearea--SexF--add50.txt", //ForSim12Oct--SexM--AverageSpeed.txt
                        "/Users/student/Documents/changearea--SexM--add50.txt"}, //ForSim30Apr--SexF--MarkovTransitionsFALSE.txt
                     (100*100), 100, // number iterations
                     "/Users/student/Documents/Sims/add50hr",
                     1000, // distance between
                     {M_PI}, // camera width
                     {55}, // camera radius
                     432, //432, 10512
                     0, //save move
                     30000,
                     0.00000001,//10000,
                     0.5 //percentage male
                     );
    
    
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
