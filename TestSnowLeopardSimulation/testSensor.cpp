//
//  testSensor.cpp
//  SnowLeopardsSimulation
//
//  Created by Elizabeth Moorcroft on 27/03/2015.
//  Copyright (c) 2015 Elizabeth Moorcroft. All rights reserved.
//

#include "catch.hpp"

#include "Sensor.h"


#include <stdio.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <random>

TEST_CASE("Sensors") {
    Sensor  sc;
    
    SECTION("interation between two lines") {
        INFO("AngleAndAngleInteraction") // Only appears on a FAIL
                std::vector <double> x_y_coord = sc.AngleAndAngleInteraction(1, 0, -1, 0 );
        CHECK(x_y_coord [0] == 0); CHECK(x_y_coord [1] == 0);
        x_y_coord = sc.AngleAndAngleInteraction(1, 0, 0.5, 0.5 );
        CHECK(x_y_coord [0] == 1); CHECK(x_y_coord [1] == 1);
        
    }
    
    //std::vector <double> Sensor::HorzAndCircInteraction(double Horz)
    SECTION("interation between circle and horizontal") {
        INFO("HorzAndAngleInteraction") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        std::vector <double> x_y_coord = sc.HorzAndCircInteraction(1);
        int check; if((x_y_coord[0] == 0 && x_y_coord[2] == 0)){check=1;} else {check=0;}
        CHECK(check == 1);
        x_y_coord = sc.HorzAndCircInteraction(0);
        if((x_y_coord[0] == -1 && x_y_coord[0] == 1) || (x_y_coord[0] == 1 && x_y_coord[2] == -1)){check=1;}else {check=0;}
        CHECK(check == 1);
        
    }
    
    //std::vector <double> Sensor::VertAndCircInteraction(double Horz)
    SECTION("interation between circle and vertical") {
        INFO("VertAndAngleInteraction") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        std::vector <double> x_y_coord = sc.VertAndCircInteraction(1);
        int check; if((x_y_coord[1] == 0 && x_y_coord[3] == 0)){check=1;} else {check=0;}
        CHECK(check == 1);
        x_y_coord = sc.VertAndCircInteraction(0);
        if((x_y_coord[1] == -1 && x_y_coord[3] == 1) || (x_y_coord[1] == 1 && x_y_coord[3] == -1)){check=1;}else {check=0;}
        CHECK(check == 1);
        
    }
    
    //AngleAndCircInteraction(double m_Angle, double c_Angle)
    SECTION("interation between circle and angle") {
        INFO("AngleAndCircInteraction") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        //AngleAndAngleInteraction(double m1_Angle, double c1_Angle, double m2_Angle, double c2_Angle)
        std::vector <double> x_y_coord = sc.AngleAndCircInteraction(1,0);
        int check; if((x_y_coord[0] == sqrt(0.5) && x_y_coord[1] == sqrt(0.5) && x_y_coord[2] == -sqrt(0.5) && x_y_coord[3] == -sqrt(0.5)) ||
                      (x_y_coord[0] == -sqrt(0.5) && x_y_coord[1] == -sqrt(0.5) && x_y_coord[2] == sqrt(0.5) && x_y_coord[3] == sqrt(0.5))){check=1;} else {check=0;}
        CHECK(check == 1);
        
    }
    //double VertAndAngleInteraction
    SECTION("interation between horizontal and angled line") {
        INFO("HorzAndAngleInteraction") // Only appears on a FAIL
        //AngleAndAngleInteraction(double m1_Angle, double c1_Angle, double m2_Angle, double c2_Angle)
        double ycoord = sc.HorzAndAngleInteraction(1, 2, 0);
        CHECK(ycoord ==0.5);
    }
    
    SECTION("interation between horizontal and angled line") {
        INFO("VertAndAngleInteraction") // Only appears on a FAIL
        double ycoord = sc.VertAndAngleInteraction(1, 2, 0);
        CHECK(ycoord ==2);
    }
    
    SECTION("distance between two points") {
        INFO("DistTwoPoints") // Only appears on a FAIL
        double dist = sc.DistTwoPoints(0, 3, 0, 4);
        CHECK(dist ==5);
    }
    
    //ngleTwoPoints(double X1, double X2, double Y1, double Y2)
    SECTION("angle between two points") {
        INFO("DistTwoPoints") // Only appears on a FAIL
        double ang = sc.AngleTwoPoints(0, 1, 0, 1);
        CHECK(ang == M_PI/4);
    }
    
    //TimeAndAngleCal(double Y, double X, double previous_y_animal, double previous_x_animal, double disttotal)
    SECTION("time and angle of capture") {
        INFO("T&AofPoints") // Only appears on a FAIL
        std::vector<double> ang = sc.TimeAndAngleCal(1, 1, 0, 0, 2*sqrt(2));
        CHECK(ang[0] == 0.5);
        CHECK(ang[1] == M_PI_4);
    }
    
    //calculate_sensoredges
    SECTION("calculate the parameters of the sensoredges") {
        INFO("calculate_sensoredges") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(0); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        CHECK(sc.get_lhs_angle() == 3*M_PI_2);
        CHECK(sc.get_rhs_angle() == M_PI_2);
        CHECK(isnan(sc.get_lhs_gradient())); CHECK(isnan(sc.get_lhs_intercept()));
        CHECK(isnan(sc.get_rhs_gradient())); CHECK(isnan(sc.get_rhs_intercept()));
        CHECK(sc.get_lhs_vert() == 0); CHECK(sc.get_lhs_horz() == 1);
        CHECK(sc.get_rhs_vert() == 0); CHECK(sc.get_rhs_horz() == 1);
    }
    
    //calculate_sensoredges
    SECTION("calculate the parameters of the sensoredges _ at angle") {
        INFO("calculate_sensoredges") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(0); sc.set_angle_halfwidth(M_PI_4);
        sc.calculate_sensoredges();
        CHECK(sc.get_lhs_angle() == 7*M_PI_4);
        CHECK(sc.get_rhs_angle() == M_PI_4);
        CHECK(sc.get_lhs_gradient() == Approx(-1.0)); CHECK(sc.get_lhs_intercept() == Approx(0));
        CHECK(sc.get_rhs_gradient() == Approx(1.0)); CHECK(sc.get_rhs_intercept() == Approx(0));
        CHECK(isnan(sc.get_lhs_vert())); CHECK(isnan(sc.get_lhs_horz()));
        CHECK(isnan(sc.get_rhs_vert())); CHECK(isnan(sc.get_rhs_horz()));
    }
    
    //calculate_sensoredges
    SECTION("calculate the parameters of the sensoredges _ vertical, direction of camera!=0") {
        INFO("calculate_sensoredges") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        CHECK(sc.get_lhs_angle() == 0);
        CHECK(sc.get_rhs_angle() == M_PI);
        CHECK(isnan(sc.get_lhs_gradient())); CHECK(isnan(sc.get_lhs_intercept()));
        CHECK(isnan(sc.get_rhs_gradient())); CHECK(isnan(sc.get_rhs_intercept()));
        CHECK(sc.get_lhs_vert() == 1); CHECK(sc.get_lhs_horz() == 0);
        CHECK(sc.get_rhs_vert() == 1); CHECK(sc.get_rhs_horz() == 0);
    }
    
    //GradientFromAngle(double angle);
    SECTION("calculting the gradient of a line from its Angle") {
        INFO("GradientFromAngle") // Only appears on a FAIL
        CHECK(sc.GradientFromAngle(M_PI/2) == Approx(0));
        CHECK(sc.GradientFromAngle(M_PI/4) == Approx(1.0));
        CHECK(sc.GradientFromAngle(3*M_PI/4) == Approx(-1.0));
    }
    
    //

    SECTION("SensorEdgeAndMovement - is there a capture when there should be") {
        INFO("SensorEdgeAndMovement- angle") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        double m_detector = sc.get_rhs_gradient();
        double g_detector = sc.get_rhs_angle();
        double c_detector = sc.get_rhs_intercept();
        double vert = sc.get_rhs_vert();
        double horz = sc.get_rhs_horz();
        sc.SensorEdgeAndMovement(2,2,//double location_x_animal, double location_y_animal,
                                -2,-2,//            double previous_x_animal, double previous_y_animal,
                                1,   //        int Individual_ID,
                                M_PI_4,    //      double move_angle,
                                 1,       //      int itnumber,
                                 1,0 ,      //      double m_animal,double c_animal,
                                 m_detector,c_detector, g_detector,  vert,  horz,     //      double m_detector, double c_detector,double g_detector, int vert, int horz,
                                 sqrt(32),Captures,1 ) ;   //      double disttotal,std::ofstream &Captures, int stepnumber)
        CHECK(sc.get_capture_current() == 1);
        CHECK(sc.get_capture_x() == 0);
        CHECK(sc.get_capture_y() == 0);
        CHECK(sc.get_capture_t() == 0.5);
        CHECK(sc.get_capture_a() == M_PI_4);
    }
    
    SECTION("SensorEdgeAndMovement - is there no capture ") {
        INFO("SensorEdgeAndMovement - no cap") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        double m_detector = sc.get_rhs_gradient();
        double g_detector = sc.get_rhs_angle();
        double c_detector = sc.get_rhs_intercept();
        double vert = sc.get_rhs_vert();
        double horz = sc.get_rhs_horz();
        sc.SensorEdgeAndMovement(20,20,//double location_x_animal, double location_y_animal,
                                 -20,-20,//            double previous_x_animal, double previous_y_animal,
                                 1,   //        int Individual_ID,
                                 M_PI_4,    //      double move_angle,
                                 1,       //      int itnumber,
                                 1,0 ,      //      double m_animal,double c_animal,
                                 m_detector,c_detector, g_detector,  vert,  horz,     //      double m_detector, double c_detector,double g_detector, int vert, int horz,
                                 sqrt(32),Captures,1 ) ;   //      double disttotal,std::ofstream &Captures, int stepnumber)
        CHECK(sc.get_capture_current() == 0);
        CHECK(isnan(sc.get_capture_x()));
        CHECK(isnan(sc.get_capture_y()));
        CHECK(isnan(sc.get_capture_t()));
        CHECK(isnan(sc.get_capture_a()));
    }
    
    SECTION("SensorEdgeAndMovement - vertical movement") {
        INFO("SensorEdgeAndMovement - vert") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        double m_detector = sc.get_rhs_gradient();
        double g_detector = sc.get_rhs_angle();
        double c_detector = sc.get_rhs_intercept();
        double vert = sc.get_rhs_vert();
        double horz = sc.get_rhs_horz();
        sc.SensorEdgeAndMovement(0,1,//double location_x_animal, double location_y_animal,
                                 0,-1,//            double previous_x_animal, double previous_y_animal,
                                 1,   //        int Individual_ID,
                                 0,    //      double move_angle,
                                 1,       //      int itnumber,
                                 NAN,NAN ,      //      double m_animal,double c_animal,
                                 m_detector,c_detector, g_detector,  vert,  horz,     //      double m_detector, double c_detector,double g_detector, int vert, int horz,
                                 2,Captures,1 ) ;   //      double disttotal,std::ofstream &Captures, int stepnumber)
        CHECK(sc.get_capture_current() == 1);
        CHECK(sc.get_capture_x() == 0);
        CHECK(sc.get_capture_y() == -1);
        CHECK(sc.get_capture_t() == 0);
        CHECK(sc.get_capture_a() == 0);
    }
    
    SECTION("SensorEdgeAndMovement - horizontal movement") {
        INFO("SensorEdgeAndMovement - horz") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        double m_detector = sc.get_lhs_gradient();
        double g_detector = sc.get_lhs_angle();
        double c_detector = sc.get_lhs_intercept();
        double vert = sc.get_lhs_vert();
        double horz = sc.get_lhs_horz();
        sc.SensorEdgeAndMovement(1,0.5,//double location_x_animal, double location_y_animal,
                                 -1,0.5,//            double previous_x_animal, double previous_y_animal,
                                 1,   //        int Individual_ID,
                                 M_PI_2,    //      double move_angle,
                                 1,       //      int itnumber,
                                 NAN,NAN ,      //      double m_animal,double c_animal,
                                 m_detector,c_detector, g_detector,  vert,  horz,     //      double m_detector, double c_detector,double g_detector, int vert, int horz,
                                 2,Captures,1 ) ;   //      double disttotal,std::ofstream &Captures, int stepnumber)
        CHECK(sc.get_capture_current() == 1);
        CHECK(sc.get_capture_x() == 0);
        CHECK(sc.get_capture_y() == 0.5);
        CHECK(sc.get_capture_t() == 0.5);
        CHECK(sc.get_capture_a() == M_PI_2);
    }
    
    //Sesnor circ
    SECTION("SensorCircAndMovement - horizontal movement") {
        INFO("SensorCircAndMovement - horz") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.SensorCircAndMovement(2,0,//double location_x_animal, double location_y_animal,
                                -2,0,// double previous_x_animal, double previous_y_animal,
                                1,// int Individual_ID,
                                M_PI_2,// double move_angle,
                                1,// int itnumber,
                                NAN,NAN,// double m_animal, double c_animal,
                                4,// double disttotal,
                                 Captures, 1);// std::ofstream &Captures, int stepnumber) ;
        int check=0;
        if((sc.get_capture_current() == 1 && sc.get_capture_x() == 1 && sc.get_capture_y() == 0 && sc.get_capture_t() == 0.75 && sc.get_capture_a() == M_PI_2) ||
           (sc.get_capture_current() == 1 && sc.get_capture_x() == -1 && sc.get_capture_y() == 0 && sc.get_capture_t() == 0.25 && sc.get_capture_a() == M_PI_2)){
            check=1;
        };
        CHECK(check == 1);
    }
    
    SECTION("SensorCircAndMovement - vertical movement") {
        INFO("SensorCircAndMovement - vert") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.SensorCircAndMovement(0,2,//double location_x_animal, double location_y_animal,
                                 0,-2,// double previous_x_animal, double previous_y_animal,
                                 1,// int Individual_ID,
                                 0,// double move_angle,
                                 1,// int itnumber,
                                 NAN,NAN,// double m_animal, double c_animal,
                                 4,// double disttotal,
                                 Captures, 1);// std::ofstream &Captures, int stepnumber) ;
        int check=0;
        if((sc.get_capture_current() == 1 && sc.get_capture_x() == 0 && sc.get_capture_y() == 1 && sc.get_capture_t() == 0.75 && sc.get_capture_a() == 0) ||
           (sc.get_capture_current() == 1 && sc.get_capture_x() == 0 && sc.get_capture_y() == -1 && sc.get_capture_t() == 0.25 && sc.get_capture_a() == 0)){
            check=1;
        };
        CHECK(check == 1);
    }
    
    SECTION("SensorCircAndMovement - angled movement") {
        INFO("SensorCircAndMovement - angle") // Only appears on a FAIL
        sc.set_radius(sqrt(2));
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.SensorCircAndMovement(2,2,//double location_x_animal, double location_y_animal,
                                 -2,-2,// double previous_x_animal, double previous_y_animal,
                                 1,// int Individual_ID,
                                 M_PI_4,// double move_angle,
                                 1,// int itnumber,
                                 1,0,// double m_animal, double c_animal,
                                 4*sqrt(2),// double disttotal,
                                 Captures, 1);// std::ofstream &Captures, int stepnumber) ;
        int check=0;
        if((sc.get_capture_current()==1 && sc.get_capture_x()==1 && sc.get_capture_y()==1 && sc.get_capture_a()==M_PI_4 && sc.get_capture_t()==Approx(0.75)) ||
           (sc.get_capture_current()==1 && sc.get_capture_x()==-1 && sc.get_capture_y()==-1 && sc.get_capture_t()==Approx(0.25) && sc.get_capture_a()==M_PI_4)
           ){
            check=1;
        };
        CHECK(check == 1);
    }
    
    SECTION("Captures inside camera range") {
        INFO("Angel to and from camera") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.set_capture_x(0);
        sc.set_capture_y(0.5);
        sc.set_capture_t(0.5);
        sc.set_capture_a(M_PI_4);
        sc.CapturesInsideCameraAngle(1,//int Individual_ID,
                                  M_PI_4,//double move_angle,
                                  1, //int itnumber,
                                  Captures,//std::ofstream &Captures,
                                  1//int stepnumber
                                     );
        CHECK(sc.get_capture_a_c2a() == 3*M_PI_2);
        CHECK(sc.get_capture_a_a2c() == 3*M_PI_4);
    }
    
    SECTION("Captures inside camera range- outside range") {
        INFO("Angel to and from camera") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.set_capture_x(-0.5);
        sc.set_capture_y(0.5);
        sc.set_capture_t(0.5);
        sc.set_capture_a(M_PI_4);
        sc.CapturesInsideCameraAngle(1,//int Individual_ID,
                                     M_PI_4,//double move_angle,
                                     1, //int itnumber,
                                     Captures,//std::ofstream &Captures,
                                     1//int stepnumber
                                     );
        CHECK(isnan(sc.get_capture_a_c2a()));
        CHECK(isnan(sc.get_capture_a_a2c()));
    }
    
    
    SECTION("Captures from (x1,y1) to (x2,y2)") {
        INFO("Captures Intersection - angle") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(2, 2, -2, -2,  1,  M_PI_4,  1, Captures,1);
        CHECK(sc.get_capture_count()==3);
    }
    
    SECTION("Captures from (x1,y1) to (x2,y2)") {
        INFO("Captures Intersection - angle2") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(2, 1, -1, -2,  1,  M_PI_4,  1, Captures,1);
        CHECK(sc.get_capture_count()==3);
    }
    
    SECTION("CapturesIntersection - vertical movement") {
        INFO("Captures Intersection - angle3") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(1, 2, 1, -2,  1,  0,  1, Captures,1);
        CHECK(sc.get_capture_count()==2);
    }
    
    SECTION("CapturesIntersection - horizontal movement") {
        INFO("Captures Intersection - angle4") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(2, 0, -2, 0,  1,  M_PI_2,  1, Captures,1);
        CHECK(sc.get_capture_count()==3);
    }
    
    SECTION("CapturesIntersection - Angle missing camera horz") {
        INFO("Captures Intersection - angle5") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(2, 3, -2, 3,  1,  M_PI_2,  1, Captures,1);
        CHECK(sc.get_capture_count()==0);
    }
    
    SECTION("CapturesIntersection - Angle missing camera vert") {
        INFO("Captures Intersection - angle5") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(3, 2, 3, -2,  1,  0,  1, Captures,1);
        CHECK(sc.get_capture_count()==0);
    }
    
    SECTION("CapturesIntersection - Angle missing camera angle") {
        INFO("Captures Intersection - angle5") // Only appears on a FAIL
        sc.set_radius(1);
        sc.set_location_x(0); sc.set_location_y(0);
        sc.set_angle_direction(M_PI_2); sc.set_angle_halfwidth(M_PI_2);
        sc.calculate_sensoredges();
        sc.set_inoutvector(1);
        sc.set_capture_count(0);
        std::ofstream CapturesNotRef;
        std::ofstream &Captures = CapturesNotRef;
        sc.CapturesIntersection(4, 8, 0, 4,  1,  0,  1, Captures,1);
        CHECK(sc.get_capture_count()==0);
    }
}
