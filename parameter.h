//
// Created by aos on 17-3-8.
//

#ifndef GPS_INS_VARIABLE_H
#define GPS_INS_VARIABLE_H
struct Acceleration{
    double x;
    double y;
    double z;
    long t;
};
struct Angle{
    double x;
    double y;
    double z;
};
struct AngularVelocity{
    double x;
    double y;
    double z;
    long t;
};
struct Position{
    double x;
    double y;
    double z;
    long t;
};
struct Velocity{
    double x;
    double y;
    double z;
    long t;
};
struct STime
{
    unsigned char ucYear;
    unsigned char ucMonth;
    unsigned char ucDay;
    unsigned char ucHour;
    unsigned char ucMinute;
    unsigned char ucSecond;
    unsigned short ucMiliSecond;
};
struct SAcc
{
    short a[3];
    short T;
};
// add Ins Vcc data struct by lc
struct SVcc
{
    short v[3];
    short T;
};
// add Ins Pos data struct by lc
struct SPos
{
    short p[3];
    short T;
};
struct SGyro
{
    short w[3];
    short T;
};
struct SAngle
{
    short Angle[3];
    short T;
};
struct SMag
{
    short h[3];
    short T;
};
struct SLonLat
{
    long lLon;
    long lLat;
};
//add quaternion data by lc
struct SQUAT
{
    short sq0;
    short sq1;
    short sq2;
    short sq3;
};
//add GPSstatus data by lc
struct SDOP
{
    short SN;
    short sPDOP;
    short sHDOP;
    short sVDOP;
};

struct SGPSV
{
    short sGPSHeight;
    short sGPSYaw;
    long lGPSVelocity;
};

struct Quaternion
{
    double q1;
    double q2;
    double q3;
    double q4;
};
#endif //GPS_INS_VARIABLE_H
