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
    std::cout << "Hello, World!\n";
    
    World sc;
    sc.MultipleIterations(100, //No cameras
                          1, //No iterations
                          "/Users/student/Documents/xx",
                          1000, //space between camera
                          {M_PI}, // camera width
                          {10}, //camera radius
                          10, //No steps
                           1, //Save Movement
                           20000, //radius
                          15000, //Max distance moved
                           2,//No states
                            {0.5,0.5}, // prob
                           {{0.5,1},{0.75,2}}, //mean
                          {{0.5,0,0,2},{0.5,0,0,2}}, ///variance
                          5 // No of animals
                          );

    
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
