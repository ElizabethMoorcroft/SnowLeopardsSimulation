//
//  Sensor.h
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 27/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#ifndef __SnowLeopardsSimulation__Sensor__
#define __SnowLeopardsSimulation__Sensor__

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

class Sensor{

private:
    int SensorID;
    
    //sensor parameters
    double location_x;
    double location_y;
    double radius;
    double angle_direction;
    double angle_halfwidth;
    
    //parameters of the sensor edges, if set up like: <)
    double lhs_angle;
    double lhs_gradient;
    double lhs_intercept;
    double lhs_vert; double lhs_horz;
    double rhs_angle;
    double rhs_gradient;
    double rhs_intercept;
    double rhs_vert; double rhs_horz;
    
    //captures
    int capture_current;
    int capture_count;
    double capture_x;
    double capture_y;
    double capture_t;
    double capture_a;
    double capture_d;
    double capture_a_a2c;
    double capture_a_c2a;
    std::vector <int> animals_in_sensor_range;
    std::vector <int> animals_just_enter_range;
    
public:
    
    void calculate_sensoredges(); //tested
    
    void set_location_x(double x){location_x=x;};
    void set_location_y(double y){location_y=y;};
    void set_radius(double r){radius=r;};
    void set_angle_direction(double a) {angle_direction=a;};
    void set_angle_halfwidth(double a) {angle_halfwidth=a;};
    
    //lhs
    double get_lhs_angle(){return (lhs_angle);};
    double get_lhs_gradient(){return(lhs_gradient);};
    double get_lhs_intercept(){return(lhs_intercept);};
    double get_lhs_vert(){return(lhs_vert);};
    double get_lhs_horz(){return(lhs_horz);};
    //rhs
    double get_rhs_angle(){return (rhs_angle);};
    double get_rhs_gradient(){return(rhs_gradient);};
    double get_rhs_intercept(){return(rhs_intercept);};
    double get_rhs_vert(){return(rhs_vert);};
    double get_rhs_horz(){return(rhs_horz);};

    //captures
    int get_capture_current(){return(capture_current);};
    double get_capture_x(){return(capture_x);};
    double get_capture_y(){return(capture_y);};
    double get_capture_t(){return(capture_t);};
    double get_capture_a(){return(capture_a);};
    double get_capture_d(){return(capture_d);};
    double get_capture_a_a2c(){return(capture_a_a2c);};
    double get_capture_a_c2a(){return(capture_a_c2a);};
    int get_capture_count(){return(capture_count);}
    void set_capture_x(double temp){capture_x = temp;};
    void set_capture_y(double temp){capture_y = temp;};
    void set_capture_t(double temp){capture_t = temp;};
    void set_capture_a(double temp){capture_a = temp;};
    void set_capture_d(double temp){capture_d= temp;};
    void set_capture_count(double temp){capture_count=temp;}
    
    void set_inoutvector(int NoAnimal){
        std::vector<int> temp(NoAnimal,0);
        animals_in_sensor_range =temp;
        animals_just_enter_range =temp;
    };

    
    
    double DistTwoPoints(double X1, double X2, double Y1, double Y2); //tested
    double AngleTwoPoints(double X1, double X2, double Y1, double Y2); //tested
    double RangeAngle(double angle); //tested
    double GradientFromAngle(double angle); //tested
    
    std::vector <double> AngleAndAngleInteraction(double m1_Angle, double c1_Angle, double m2_Angle, double c2_Angle); //tested
    double VertAndAngleInteraction(double Vert, double m_Angle, double c_Angle); //tested
    double HorzAndAngleInteraction(double Horz, double m_Angle, double c_Angle); //tested
    std::vector <double> HorzAndCircInteraction(double Horz); //tested
    std::vector <double> VertAndCircInteraction(double Vert); //tested
    std::vector <double> AngleAndCircInteraction(double m_Angle, double c_Angle); //tested
    std::vector <double> TimeAndAngleCal(double Y, double X, double previous_y_animal, double previous_x_animal, double disttotal); //tested
    
    
    void SensorEdgeAndMovement(double location_x_animal, double location_y_animal,double previous_x_animal, double previous_y_animal,int Individual_ID,  double move_angle, int itnumber, double m_animal,double c_animal, double m_detector, double c_detector,double g_detector,double vert, double horz, double disttotal,std::ofstream & Captures, int stepnumber); //tested
    void SensorCircAndMovement(double location_x_animal, double location_y_animal,  double previous_x_animal, double previous_y_animal, int Individual_ID, double move_angle, int itnumber,double m_animal, double c_animal, double disttotal, std::ofstream &Captures, int stepnumber); //tested
    void CapturesInsideCameraAngle(int Individual_ID, double move_angle,int itnumber, std::ofstream &Captures,int stepnumber); //tested
    void CapturesIntersection(double location_x_animal, double location_y_animal, double previous_x_animal, double previous_y_animal, int Individual_ID,   double move_angle, int itnumber, std::ofstream &Captures, int stepnumber);
    
    void UpdateCaptures(int Animal_ID, int itnumber, std::ofstream &Captures, int stepnumber);
    void animal_in_out_range(int animal_id, double time, int stepnumber);
    
};

#endif /* defined(__SnowLeopardsSimulation__Sensor__) */
