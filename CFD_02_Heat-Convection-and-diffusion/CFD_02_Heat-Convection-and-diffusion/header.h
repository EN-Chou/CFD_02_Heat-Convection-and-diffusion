#pragma once
#include <mpi.h>
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

//Problem Description
#define PI 3.1415926
#define N 161
#define Pe numeric_limits<double>::infinity()
#define scheme 3 //0:Upwind, 1:2nd Up, 2: CDS, 3: QUICK, 4: minmod
#define SOR 1
#define BC_exchange_rate 1
#define tol pow(10, -9)


void write(double* a, int x, int y);
void printinfo();
void phi_generate(int);


















