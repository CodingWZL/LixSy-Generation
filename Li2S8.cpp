#include<stdio.h>
#include<string.h>
#include<iostream>
#include<stdlib.h>
#include<sys/time.h>
#include<math.h>
#include<unistd.h>
using namespace std;

#define pi acos(-1) // the value of pi

double center_x = 0, center_y = 0, center_z = 0;  // the center of LixSy
double low_x = 0, low_y = 0, low_z = 0;   // the xyz responding to the lowest value of z


/************
 * the structure of atom, including xyz, radial, chemical special
 * the first 10 atoms belongs to Li2S8
 * the last 48 atoms belongs to MoSe2
************/
struct E{
    double x,y,z;
    double radial;
    string special;
}atoms[60];


/************
 * Product the random number, including the rotation angles and random points of xyz
************/
double getRandData(double low, double high)
{     
      struct timeval tv;
      gettimeofday(&tv, NULL);
      srand(tv.tv_sec + tv.tv_usec + getpid());
      double m=rand()/(double)(RAND_MAX) * (high - low) + low;
      return m;
}


/************
 * Read the xyz in the POSCAR.Li2S8 and POSCAR.MoSe2
************/
void readxyz(int atom)
{
    int k,j=8;
    FILE *init;
    double tmp1, tmp2, tmp3;
    char StrLine[1024];

    // read the Li2S6
    init = fopen("POSCAR-Li2S8", "r");
    while(j--)
    {
        fgets(StrLine, 1024, init);
    }
    for(k = 1; k <= atom; k++)
    {
        fscanf(init, "%lf%lf%lf", &atoms[k].x, &atoms[k].y, &atoms[k].z);
        if(k<=2)
        {
            atoms[k].special = "Li";
            atoms[k].radial = 1.57;
        }
        else
        {
            atoms[k].special = "S";
            atoms[k].radial = 1.04;
        }
        
        // Convert fractional coordinates to direct coordinates
        tmp1 = atoms[k].x;
        tmp2 = atoms[k].y;
        tmp3 = atoms[k].z;
        atoms[k].x = tmp1 * 15.75777018 - tmp2 * 0.03238349 + tmp3 * 0.67693482;
        atoms[k].y = tmp2 * 15.79031032 - tmp1 * 0.02240748 - tmp3 * 0.25966241;
        atoms[k].z = tmp3 * 13.19026268 + tmp1 * 0.674226290 - tmp2 * 0.21974230;
    }
    fclose(init);

    // read the MoSe2
    init = fopen("POSCAR-MoSe2", "r");
    j=8;
    while(j--)
    {
        fgets(StrLine, 1024, init);
    }
    for(k = atom+1; k <= atom+48; k++)
    {
        fscanf(init, "%lf%lf%lf", &atoms[k].x, &atoms[k].y, &atoms[k].z);
        if(k<=22)
        {
            atoms[k].special = "Mo";
            atoms[k].radial = 1.40;
        }
        else
        {
            atoms[k].special = "Se";
            atoms[k].radial = 1.04;
        }
       
        // Convert fractional coordinates to direct coordinates
        atoms[k].x = atoms[k].x * 13.2841 - atoms[k].y * 6.64205;
        atoms[k].y = atoms[k].y * 11.504367;
        atoms[k].z = atoms[k].z * 23.337;
    }
    fclose(init);

}

// Calculate the center coordination of Li2S8
void center(int atom)
{
    int k;
    for(k = 1;k <= atom;k++)
    {
        center_x = center_x + atoms[k].x;
        center_y = center_y + atoms[k].y;
        center_z = center_z + atoms[k].z;
    }
    center_x = center_x/atom;
    center_y = center_y/atom;
    center_z = center_z/atom;
}

/**************
 * The rotation matrix
**************/
void rotation(int atom, double alpha, double beta, double gamma)
{
    double tmp;
    int k;
    low_z = 100;
    for(k = 1;k <= atom;k++)
    {
        // Fisrtly, shift the center of rotation to the origin
        atoms[k].x = atoms[k].x - center_x;
        atoms[k].y = atoms[k].y - center_y;
        atoms[k].z = atoms[k].z - center_z;

        //X rotation
        tmp = atoms[k].y;
        atoms[k].x = atoms[k].x;
        atoms[k].y = atoms[k].y*cos(alpha*pi)+atoms[k].z*sin(alpha*pi);
        atoms[k].z = (-1.0)*tmp*sin(alpha*pi)+atoms[k].z*cos(alpha*pi);

        //Y rotation
        tmp = atoms[k].x;
        atoms[k].x = atoms[k].x*cos(beta*pi)-atoms[k].z*sin(beta*pi);
        atoms[k].y = atoms[k].y;
        atoms[k].z = tmp*sin(beta*pi)+atoms[k].z*cos(beta*pi);

        //Z rotation
        tmp = atoms[k].x;
        atoms[k].x = atoms[k].x*cos(gamma*pi)+atoms[k].y*sin(gamma*pi);
        atoms[k].y = (-1.0)*tmp*sin(gamma*pi)+atoms[k].y*cos(gamma*pi);
        atoms[k].z = atoms[k].z;
        
        // Shift it back
        atoms[k].x = atoms[k].x + center_x;
        atoms[k].y = atoms[k].y + center_y;
        atoms[k].z = atoms[k].z + center_z;
    }

    // Calculate the xyz which has the lowest value of z
    for(k = 1;k <= atom;k++)
    {
        if(low_z >= atoms[k].z)
        {
            low_x = atoms[k].x;
            low_y = atoms[k].y;
            low_z = atoms[k].z;
        }
    }
}

/***************
 * Write the resulting coordinates to POSCAR
***************/
void writexyz(int atom, double translation_x, double translation_y, double translation_z)
{
    FILE *prod;
    int a,k;
    string base_path = "POSCAR";
    prod = fopen(base_path.c_str(), "a+");
    fprintf(prod, "CONTCAR\\(0\\0\\1)\n1.0\n");
    fprintf(prod, "      13.2840995788       0.0000000000       0.0000000000\n");
    fprintf(prod, "      -6.6420497893      11.5043677016       0.0000000000\n");
    fprintf(prod, "       0.0000000000       0.0000000000      23.3369998932\n");
    fprintf(prod, "     Mo   Se   Li   S\n");
    fprintf(prod, "     16   32   2   %d\n", atom-2);
    fprintf(prod, "Cartesian\n");

    for(k = atom+1; k <= atom+48;k++)
    {
        fprintf(prod, "      %.8f     %.8f     %.8f\n", atoms[k].x, atoms[k].y, atoms[k].z);
    }
    for(k = 1; k <= atom; k++)
    {
        fprintf(prod, "      %.8f     %.8f     %.8f\n", atoms[k].x+translation_x, atoms[k].y+translation_y, atoms[k].z+translation_z);
       // printf("      %.8f     %.8f     %.8f\n", atoms[k].x+translation_x, atoms[k].y+translation_y, atoms[k].z+translation_z);
    }
    fclose(prod);
}

/**************
 * Main function
**************/
int main()
{
    int atom = 10;   // No. of atoms of Li2S8
    double alpha, beta, gamma;   // rotation angle
    double rand_x, rand_y, rand_z;  // random point
    double translation_x, translation_y, translation_z;  //translation value

    
    alpha = getRandData(0, 2.0);   // the rotation angle
    beta = getRandData(0, 2.0);
    gamma = getRandData(0, 2.0);

    readxyz(atom);  // read Li2S8 and MoSe2
    center(atom);   // calculate the center of Li2S8

    rotation(atom, alpha, beta, gamma);  // rotate the Li2S8

    // random points above the MoSe2
    rand_x = getRandData(0.0, 1.0);
    rand_y = getRandData(0.0, 1.0);
    rand_z = getRandData(0.237098, 0.262808);   //2.2+3.33317----2.8+3.33317(5.53317----6.13317)

    // direct coordinates of random point
    rand_x = rand_x * 13.2841 - rand_y * 6.64205;
    rand_y = rand_y * 11.504367;
    rand_z = rand_z * 23.337;
    //printf("%f %f %f %f %f %f %f\n", atoms[1].z, atoms[2].z, atoms[3].z, atoms[4].z, atoms[5].z, atoms[6].z, low_z);
    // shift value
    translation_x = rand_x - low_x;
    translation_y = rand_y - low_y;
    translation_z = rand_z - low_z;
    // product POSCAR
    writexyz(atom, translation_x, translation_y, translation_z);

    return 0;
}

