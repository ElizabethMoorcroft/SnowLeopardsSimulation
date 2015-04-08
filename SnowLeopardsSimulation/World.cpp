//
//  World.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 08/04/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

//Header files
#include "Animal.h"
#include "Sensor.h"
#include "World.h"



std::vector<Sensor*>  World::trappinggrid(int NoSensors, double cam_interval, std::ofstream &Sensors, std::vector<double> SensorWidth, std::vector<double> SensorRadius){
    
    int no_cam_locations_y = sqrt(NoSensors); int no_cam_locations_x = sqrt(NoSensors);

    std::vector<Sensor*> AllSensors(pow(no_cam_locations_y,2));
    
    int sensorcount=0;
    int LengthSW = SensorWidth.size();
    int LengthSR = SensorRadius.size();
    
    
    for(int grid_number_x = 0; grid_number_x < no_cam_locations_x; grid_number_x ++){
        double x_location = grid_number_x*cam_interval;
        
        for(int grid_number_y = 0; grid_number_y < no_cam_locations_y; grid_number_y ++){
            double y_location = grid_number_y*cam_interval;
            
            for(int sensor1 =0; sensor1<LengthSW; sensor1 ++){
                for(int sensor=0; sensor<LengthSR; sensor++){
                    //std::cout<< "sensor: " << sensorcount << "/"<< NoSensors <<std::endl;
                    AllSensors[sensorcount] =new Sensor(sensorcount, x_location, y_location, SensorRadius[sensor],0, SensorWidth[sensor1]);
                    //Saves the locations and the angle of the Sensor
                    Sensors << sensorcount << //1st column
                        "," << x_location << //2nd column
                        "," << y_location << //...
                        "," << 0 << //...
                        "," << SensorWidth[sensor1] << //4th column
                        "," << SensorRadius[sensor] << //4th column
                        "," << "Grid" <<
                        "," << grid_number_x <<
                        "," << grid_number_y <<
                        "\n";
                    sensorcount +=1;
                };
            };
        };
    };
    return AllSensors;
};


void World::AnimalMovement(int id, double randomstart,
                           int Iteration, int SaveMovement,
                           double CentreHome_r,double MaximumDistance,
                           double no_of_move_states,std::vector<double> probability,
                           std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
                           int NoSteps, std::ofstream &Movement ,  std::vector<Sensor*> AllSensors , std::ofstream &Captures,
                           double cam_interval,int  NoSensors){
    srand(randomstart);
    double seed1, seed2, seed3,seed4,seed5, seed6, temp;
    for(int j=0; j<251; j++){
        if(j==50){seed1=double(rand());}
        else if(j==100){seed2=double(rand());}
        else if(j==150){seed3=double(rand());}
        else if(j==200){seed4=double(rand());}
        else if(j==250){seed5=double(rand());}
        else {temp=double(rand());};
    };
    
    //sets random location
    double xlocation = StartLocation(NoSensors, cam_interval,seed1, CentreHome_r);
    double ylocation = StartLocation(NoSensors, cam_interval,seed2, CentreHome_r);
    
    // To choose a start angle, sets up a random number class
    // Uses a radom number from stream RandomNumberStreamAnimalAngle for a random seed
    srand(seed3);double CurrentAngleTemp = ((double) rand()/RAND_MAX)*2*M_PI;
    
    // New animal given start locations - at the centre of the home range
    //  Inputs are: ID & Starts location (x,y) &  Initial angle
    Animal Animal1(id,Iteration, SaveMovement,
                   CentreHome_r,MaximumDistance,
                   no_of_move_states,probability,mean_vector,variance_vector,
                   xlocation, ylocation, CurrentAngleTemp, seed4);
    
    //Sets seed for a random number
    //Random number stream for the movemnet of the animal
    srand(seed5);
    for(int j=0; j<NoSteps; j++){
        seed6 =double(rand());
        for(int extra=0; extra<999; extra++){temp=double(rand());};
        Animal1.UpdateLocation(seed6, Movement,  AllSensors, Captures);
    }; //End of j loop for Steps

}

double World::StartLocation(int NoSensors, double cam_interval, double seed, double CentreHome_r){
    
    int no_locations = sqrt(NoSensors);
    double width = no_locations*cam_interval;
    
    srand(seed); double output= ((double) rand()/RAND_MAX)*(width + CentreHome_r*2)-CentreHome_r;
    
    return(output);
};




std::string World::make_filename( const std::string& directory ,const std::string& basename){
    std::ostringstream result;
    result << directory <<basename;
    return result.str();
};



void World::oneiteration(
                double cam_interval,
                int NoSteps,
                int Iteration, int SaveMovement,
                double CentreHome_r,double MaximumDistance,
                double no_of_move_states,std::vector<double> probability,
                std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector,
                double seed,
                std::ofstream &Captures, std::ofstream &Movement,std::vector<Sensor*> AllSensors, int NoAnimal){
    
    double seed1, temp;
    
    int lengthsensors = AllSensors.size();
    //run animals
    for(int i=0; i<NoAnimal; i++){
        
        srand(seed);
        for(int j=0; j<NoAnimal; j++){
            seed1=double(rand());
            for(int k=0; k<NoAnimal;k++){temp=double(rand());}
            seed= double(rand());
        };
        
        std::cout<<seed1<<std::endl;
        AnimalMovement(i, seed1, Iteration, SaveMovement,
                        CentreHome_r,MaximumDistance,no_of_move_states,probability,mean_vector,variance_vector,
                       NoSteps, Movement ,  AllSensors , Captures,cam_interval, lengthsensors);
    }

};

void World::MultipleIterations(int NoSensors, int NoInterations, std::string savevalue,
                                 double cam_interval, std::vector<double> SensorWidth, std::vector<double> SensorRadius,
                                 int NoSteps,
                                 int SaveMovement,
                                 double CentreHome_r,double MaximumDistance,
                                 double no_of_move_states,std::vector<double> probability,
                                 std::vector<std::vector<double>> mean_vector,std::vector<std::vector<double>> variance_vector, int NoAnimals 
                                 ){
    
    //Creates file for Sensors (CSV file) and writes in the header
    std::ofstream Sensors;
    Sensors.open(make_filename(savevalue, ",Sensors.csv" ).c_str());
    Sensors << "sensor.ID" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "centre.angle" <<
    "," << "half.width.angle" <<
    "," << "radius" <<
    "," << "placement" <<
    "," << "grid.row" <<
    "," << "grid.column" <<
    "\n";
    
    //Creates file for Captures (CSV file) and writes in the header
    std::ofstream CapturesNotRef;
    std::ofstream &Captures = CapturesNotRef;
    Captures.open(make_filename(savevalue, ",Captures.csv" ).c_str());
    Captures << "animal.number" <<
    "," << "step.number" <<
    "," << "sensor.ID" <<
    "," << "iteration.number" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "time" <<
    "," << "angle.from.camera.to.animal" <<
    "," << "angle.from.animal.to.camera" <<
    "," << "distance.to.camera" <<
    "," << "in.range"
    "," << "just.entered"
    "\n";
    
    //Creates file for Movement (CSV file) and writes in the header
    std::ofstream MovementNotRef;
    std::ofstream &Movement = MovementNotRef;
    Movement.open(make_filename(savevalue, ",Movement.csv" ).c_str());
    Movement << "animal.number" <<
    "," << "step.number" <<
    "," << "x.location" <<
    "," << "y.location" <<
    "," << "movement.angle" <<
    "," << "total.distance" <<
    "," << "current.state" <<
    "," << "edge.interaction" <<
    "," << "distance.moved" <<
    "," << "iteration.number" <<
    "," << "current.state" <<
    "\n";
    
    
    
    //set up grid
    std::vector<Sensor*>  AllSensors = trappinggrid(NoSensors, cam_interval, Sensors,  SensorWidth, SensorRadius);
    
 
    double seed, temp;
    for(int iteration =0; iteration<NoInterations;iteration++){
        
        srand(iteration);
        for(int j=0; j<NoInterations; j++){
            seed=double(rand());
        };
        
        oneiteration(cam_interval,
                     NoSteps,
                     iteration, SaveMovement,
                     CentreHome_r,MaximumDistance,
                     no_of_move_states, probability,
                     mean_vector, variance_vector,
                     seed,
                     Captures, Movement, AllSensors, NoAnimals);
      
    };

};

