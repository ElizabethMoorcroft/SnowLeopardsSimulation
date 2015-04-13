//
//  Animal.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 25/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include "Animal.h"
#include "Sensor.h"
#include "CheckingSensors.h"
#include <random>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>

Animal::Animal(){};
Animal::Animal(int AnimalId,int iteration,int savemovement,
               // movement setups
               double CentreHome_r,double MaximumDistance,
               double no_of_move_states,std::vector<double> probability,std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
               std::vector<std::vector<double>> transitions,
               double x, double y, double a, double seed,
               std::ofstream &Movement, std::vector<Sensor*> AllSensors , std::ofstream &Captures, int sex
               ){
    set_AnimalId(AnimalId);
    set_MaximumDistance(MaximumDistance);
    
    set_Current_x(x);
    set_Current_y(y);
    
    
    set_no_of_move_states(no_of_move_states);
    set_probability(probability);
    set_mean_vector(mean_vector);
    set_variance_vector(variance_vector);
    transitions_vector = transitions;
    
    set_CentreHome_x(x);
    set_CentreHome_y(y);
    set_CentreHome_r(CentreHome_r);
    
    int s = PickState (seed, no_of_move_states,  probability);
    set_Current_state(s);
    
    step_number = 0;
    Total_distance = 0;
    Sex = sex;
    SaveMovement = savemovement;
    Iteration = iteration;
    //std::cout<<get_Current_x()<<std::endl;
    LocationVector(x, y, 0,1, Movement, AllSensors , Captures);

};

void Animal::LocationVector(double previous_x, double previous_y, int LeaveEntreCode, int End, std::ofstream &Movement, std::vector<Sensor*> AllSensors , std::ofstream &Captures){
    
    if(SaveMovement==1 || step_number==0 || step_number == 431){
        //std::cout<< "movement being saved" <<std::endl;
        
        
        Movement<< AnimalId << //1st column, row "stepcounter"
            "," << step_number << //2nd column, row "stepcounter"
            "," << Current_x << //...
            "," << Current_y << //...
            "," << Current_angle << //...
            "," << Total_distance << //...
            "," << Current_state << //...
            "," << LeaveEntreCode << //8th column, row "stepcounter"
            "," << Current_distance << //8th column, row "stepcounter"
            "," << Iteration <<                  // itertaion number
            "," << Sex <<
            "\n";
    };
    int NoSensors = AllSensors.size();
    CheckingSensor(Current_x,Current_y, Current_angle, previous_x, previous_y,
                   AllSensors,
                   Movement,
                   Captures,
                   step_number, AnimalId, Iteration,
                   NoSensors);

};

void  Animal::UpdateLocation (double seed, std::ofstream &Movement, std::vector<Sensor*> AllSensors , std::ofstream &Captures){
    
    double previousx = Current_x;
    double previousy = Current_y;
    
    // This starts a stream of random numbers used twice
    //  -> To start a new stream of RandNum for the update of movement
    //  -> To start a new stream of RandNum for the probability of changing states in 2 state corr walk (MoveType==2)
    std::vector<double> seed_stream=ThreeSeeds(seed);
    
    //std::cout<<"seed " <<seed <<" seed_stream " <<seed_stream[0]<<std::endl;
    NewLocation(seed_stream[0]);
    
    LocationVector(previousx, previousy, 0,1, Movement, AllSensors , Captures);

    
    //Increases step number by one once the animal finishes moving in the environment
    step_number =step_number+1;
    
};//End of update location



// For cluster analyis
void Animal::NewLocation (double seed){
    
    //std::cout<<"seed "<< seed <<std::endl;
    std::vector<double> seed_stream=ThreeSeeds(seed);
    //Pick state
    int new_state = PickState(seed_stream[0], no_of_move_states, transitions_vector[Current_state]);
    std::vector<double> mean = mean_vector[new_state];
    std::vector<double> variance = variance_vector[new_state];
    std::vector<double> newlocations(2);
    
    //std::cout<<" new_state "<< new_state <<std::endl;
    
    int count =0;
    newlocations[0] = CentreHome_r+CentreHome_x; newlocations[1] = CentreHome_r+CentreHome_y;
    while(sqrt(pow(newlocations[0]-CentreHome_x,2)+pow(newlocations[1]-CentreHome_y,2))>CentreHome_r && count <(10^4)){
        //std::cout<<"seed_stream[2] "<< seed_stream[2]<<std::endl;
        seed_stream=ThreeSeeds(seed_stream[2]);
        newlocations = NewLocationFromMult(seed_stream[1], mean, variance);
        //std::cout<<" location (" <<newlocations[0]<<","<<newlocations[1]<<")"<<std::endl;
        count +=1;
        
    }
    if(count>=(10^4) && count<=(10^4)*2){
        
        while(sqrt(pow(newlocations[0]-CentreHome_x,2)+pow(newlocations[1]-CentreHome_y,2))>CentreHome_r){
            
            seed_stream=ThreeSeeds(seed_stream[2]);
            new_state = PickState(seed_stream[0], no_of_move_states, probability);
            mean = mean_vector[new_state];
            variance = variance_vector[new_state];
            newlocations = NewLocationFromMult(seed_stream[1], mean, variance);

        }
    } else{
        while(sqrt(pow(newlocations[0]-CentreHome_x,2)+pow(newlocations[1]-CentreHome_y,2))>CentreHome_r){
            
            seed_stream=ThreeSeeds(seed_stream[2]);
            newlocations = NewLocationFromMult(seed_stream[1], {5,M_PI}, {2,0,0,2} ); //these are test inputs
            
        }
    
    }
    
    Current_state = new_state;
    Current_x = newlocations[0];
    Current_y = newlocations[1];
    Current_distance = newlocations[2];
    Current_angle = newlocations[3];
    Total_distance += Current_distance;
};

std::vector<double> Animal::NewLocationFromMult(double seed, std::vector<double> mean, std::vector<double> variance){
    std::vector<double> values = MultivariteNorm(seed,mean, variance);
    values[0] = exp(values[0]);
    
    double NextDist = values[0];
    double NextAngle = values[1]+Current_angle;
    NextAngle = RangeAngle(NextAngle); //Recalculates so that the angle between 0 and 360
    std::vector<double>  newlocation = {CalNext_X(NextDist, NextAngle),CalNext_Y(NextDist, NextAngle),NextDist,NextAngle};
    return newlocation;
}

//Calculate multivarite normal
std::vector<double> Animal::MultivariteNorm (double seed, std::vector<double> mean,std::vector<double> variance){
    
    //Multivaraite normal:
    //    X = MEAN + NORMAL%*%CHOLSIGMA
    //Where:
    //      X has the dimension x
    //      NORMAL is a row vector length N
    //      CHOLSIGMA is the Cholesky transformation of the sigma matrix
    
    srand(seed);
    std::vector<double> seed_stream=ThreeSeeds(seed);
    double seed_start = seed_stream[0];
    
    // Multivariate Normal
    std::vector<double>  values = {-1,-1};
    std::vector<double>  tempvector(2);
    
   // int count=1;
    while(values[0]<=0 || values[1]<=0 || values[0]> log(MaximumDistance) || values[1] > 2*M_PI){
       // std::cout<<count<<std::endl; count+=1;
        std::vector<double> seed_loop=ThreeSeeds(seed_start);
        // Generate two random normal values
        std::vector<double>  NORMAL(2);
        for(int n=0; n<2;n++){NORMAL[n] = Normal(seed_loop[n], 0, 1);}
        tempvector = MatrixMultip(NORMAL, variance );
        
       // std::cout<<"tempvector " << tempvector[0] <<", "<< tempvector[1] <<std::endl; count+=1;
        for(int n=0; n<2;n++){values[n] = mean[n] + tempvector[n];}
        values[1] = RangeAngle(values[1]);
        seed_start = seed_loop[2];
        
    } //end of while loop - will exit when both values>0
    
    // Makes
    //if(AtoBUnif (seed_stream[2],-1,1)<0){values[1]=values[1]*-1;}
    //values[0] = exp(values[0]);
    
    return values;
};

double Animal::Normal (double seed, double Mean, double SD){
    
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution(Mean,SD);
    double number = distribution(generator);
    return number;
};

std::vector<double> Animal::MatrixMultip (std::vector<double> RANDOM, std::vector<double> variance) {
    //where the matrix:
    // ( a b )
    // ( c d )
    //is in vector form[a,c,b,d]
    //adapted from: www.code.wikia.com/wiki/Matrix_multiplication
    std::vector<double> C = {0,0};
    for (int col = 0; col < 2; ++col){
        double sum = 0;
        for (int inner = 0; inner < 2; ++inner){
            sum += RANDOM[inner] * variance[inner+col*2];
        }
        C[col] = sum;
    }
    
    return C;
};

int Animal::PickState (double Seed,  double no_of_move_states, std::vector<double> probability){
    
    //Sets seed for a random number
    srand(Seed);
    
    //double no_of_move_states=parameter_values[0][0][0];
    //std::vector<double> probability=parameter_values[3][0];
    
    // Uniform 0 - 1 random number
    double v1 = ((double) rand()/RAND_MAX);
    
    int state = -1;
    double prob = 0;
    int n=0;
    //
    while( state<0 ){
        if(n==no_of_move_states){std::cout<< "ERROR - STATE NOT FOUND"<< std::endl; exit(EXIT_FAILURE);}
        prob +=  probability[n];
        if(v1 < prob) {state = n;};
        n+=1;
        //std::cout<< "n " <<n << " no_of_move_states " <<no_of_move_states << " v1 " <<v1 <<std::endl;
    }
    
    //returns value
    return state;
}

std::vector<double> Animal::ThreeSeeds(double seed){
    srand(seed);
    double temp;
    std::vector<double> output(3);
    for(int i=0; i<(601); i++){
        if(i==200){output[0] =double (rand());}
        else if(i==400){output[1] =double (rand());}
        else if(i==600){output[2] =double (rand());}
        else {temp =double (rand());}
    };
    return output;
}
double Animal::CalNext_X(double distance, double angle){double newx = Current_x + distance*sin(angle); return newx;};
double Animal::CalNext_Y(double distance, double angle){double newy = Current_y + distance*cos(angle); return newy;};
double Animal::RangeAngle(double angle){
    while(angle<0){angle += 2*M_PI;}
    while(angle>=2*M_PI){angle -= 2*M_PI;};
    return angle;
};
