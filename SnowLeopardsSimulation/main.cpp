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

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    Animal RandomNumber1;
    int state = RandomNumber1.PickState(1, 2, {0,1});
    RandomNumber1.set_CentreHome_x(0);RandomNumber1.set_CentreHome_y(0);RandomNumber1.set_CentreHome_r(10000000);
    RandomNumber1.set_Current_x(0);RandomNumber1.set_Current_y(0);RandomNumber1.set_Current_state(0);
    RandomNumber1.set_MaximumDistance(15000);
    for(int j=0; j<20; j++){
        //std::cout   << "j: "<< j <<std::endl;
        RandomNumber1.NewLocation(j, 2, {0.5,0.5},{{1,0},{10,M_PI}},{{0.1,0,0,0.01},{1,0,0,0.1}});
        std::cout   << "j: "<< j
                    <<", state-"<<RandomNumber1.get_Current_state()
                    <<" location (" <<RandomNumber1.get_Current_x()<<","<<RandomNumber1.get_Current_y()<<")"
                    <<std::endl;
        
    }
    
    //std::cout << "state" <<state <<std::endl;
    //std::cout << "returnedvector (" <<returnedvector[0]<<","<<returnedvector[1]<<")"<<std::endl;
    return 0;
}
