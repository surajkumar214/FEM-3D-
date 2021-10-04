#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#include<stdio.h>
//to evaluate the determinants
float Determinants(float [4][4], int i, int j);
//matrix multiplication
void Matrix_Multiply(float [6][6] , float[6][12] ,float [6][12] );
// transpose and matrix multiplication
void Trans_Matrix_Multiply(float  [6][12], float [6][12], float [12][12]);
//evaluation of B matrix
void B_matrix(float [4][4],float [6][12]);


#define N 4   // no of nodes
#define E 1  // no of elements



std::vector<float> positions;
std::vector<unsigned int> indices;
//////////////////////////////////////////////////////parser
void parse_file(const std::string &filename)
{


    std::ifstream fin(filename, std::ios::in);
    std::string lineBuffer;
    bool vert_parse = false, elem_parse = false,start = false;
    // float myvec[3];
    // int myind[4];
    while (std::getline(fin, lineBuffer))
    {
        std::stringstream ss(lineBuffer);
        std::string cmd;
        ss >> cmd;

 if (vert_parse )
        {
            // glm::vec3 vertex;
            // int dim = 0;
            float dummy;
            //ss >> dummy;
            float x, y, z;
            ss >> x;
            ss >> y;
            ss >> z;

            std::cout << x << " " << y << " " << z << std::endl;

            positions.push_back(x);
            positions.push_back(y);
            positions.push_back(z);
        }
        else if (elem_parse)
        {
            int dummy;
            //ss >> dummy;
            unsigned int x, y, z, w;
            ss >> x;
            ss >> y;
            ss >> z;
            ss >> w;
            std::cout << x << " " << y << " " << z << " " << w << std::endl;

            indices.push_back(x);
            indices.push_back(y);
            indices.push_back(z);
            indices.push_back(w);
        }

        if (cmd == "*VERTICES")
        {
            vert_parse = true;
            std::getline(fin, lineBuffer);

        }
        if (cmd == "*ELEMENTS")
        {
            elem_parse = true;
            vert_parse = false;
            std::getline(fin, lineBuffer);
            std::getline(fin, lineBuffer);

        }

    }
    // std::cout << positions[0] << " " << positions[1] << " " << positions[2] << std::endl;
    // std::cout << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3] << std::endl;
}

int main()
{
    parse_file("live_c_unity.veg");





   for(int i=0;i<E;i++)
     {

       int posit = 4*i;

	int v1 = indices[posit]; int v2 = indices[posit+1]; int v3 = indices[posit+2]; int v4 = indices[posit+3];
 	int tets [4]; tets[0] = v1;tets[1] =v2; tets[2] = v3; tets[3] = v4;





       float Coordinates[4][4];
Coordinates[0][0] = 1; Coordinates[0][1] = positions[3*v1];Coordinates[0][2] = positions[3*v1+1];; Coordinates[0][3]=positions[3*v1+2];
Coordinates[1][0] = 1; Coordinates[1][1] = positions[3*v2];Coordinates[1][2] = positions[3*v2+1];; Coordinates[1][3]=positions[3*v2+2];
Coordinates[2][0] = 1; Coordinates[2][1] = positions[3*v3];Coordinates[2][2] = positions[3*v3+1];; Coordinates[2][3]=positions[3*v3+2];
Coordinates[3][0] = 1; Coordinates[3][1] = positions[3*v4];Coordinates[3][2] = positions[3*v4+1];; Coordinates[3][3]=positions[3*v4+2];



       int tetid = i;

////////////////////////////////////////////////////////////local stiffness
float Shape_function[4][4];
static float B [6][12];


static float K[12][12],k1[6][12];

static float C[6][6];
float Y = 1, v = 0;
float l = (Y*v)/((1+v)*(1-2*v)),u = Y/((1+v));

C[0][0] =l+u;
C[0][1] =l;
C[0][2] = l;

C[1][1] = l+u;
C[1][0] = l;
C[1][2] = l;



C[2][2]= l+u;
C[2][0] = l;
C[2][1] = l;


C[3][3] = u/2;C[4][4]=u/2;C[5][5]=u/2;


//we get 6 times the volume, 60 floating point operations
float Volume6 = Determinants(Coordinates, 0,0)-Determinants(Coordinates, 1,0)+Determinants(Coordinates, 2,0)-Determinants(Coordinates, 3,0);
printf("%f",Volume6);
// evaluation of shape function, 224 floating point operations
for(int i =0;i<4;i++)
for(int j=0;j<4;j++)
{
float det = Determinants(Coordinates, i,j);
if((i+j)%2==0)
Shape_function[i][j] = det;
else
Shape_function[i][j] = -1*det;

}




// evaluation of the constant B matrix 0 floating point operations
B_matrix(Shape_function , B);
printf("\n  \n");


//432*2 floating point operations
Matrix_Multiply( C ,B,k1);
printf("\n");


//864*2 floating point operations
Trans_Matrix_Multiply(B,k1, K);

//144 floating operations
for(int i =0;i<12;i++)
for(int j=0;j<12;j++)
{
K[i][j] = (K[i][j])/(6*Volume6);
}



for(int i=0;i<N*3;i++)
     {
          for(int j=0;j<N*3;j++)
          {
               if (K[i][j]!=0)
               {
               printf("%d  %d  %f \n",i,j,K[i][j]);

               }

          }

     }


/////////////////////////////////////////////////////////// global  assembling







}






return 0;
}

/////////////////////////////////////////////
//////////////////////////////////////////////////////calculation of B matrix

void  B_matrix(float shape[4][4],float matrix[6][12])
{

 matrix[0][0] = shape[0][1]; matrix[0][3] = shape[1][1];  matrix[0][6] = shape[2][1];matrix[0][9] = shape[3][1];
 matrix[1][1] = shape[0][2]; matrix[1][4] = shape[1][2];  matrix[1][7] = shape[2][2];matrix[1][10] = shape[3][2];
 matrix[2][2] = shape[0][3]; matrix[2][5] = shape[1][3];  matrix[2][8] = shape[2][3];matrix[2][11] = shape[3][3];

 matrix[3][0] = shape[0][2]; matrix[3][1] = shape[0][1]; matrix[3][3] = shape[1][2]; matrix[3][4] = shape[1][1];
 matrix[3][6] = shape[2][2]; matrix[3][7] = shape[2][1]; matrix[3][9] = shape[3][2]; matrix[3][10] = shape[3][1];

 matrix[4][1] = shape[0][3]; matrix[4][2] = shape[0][2]; matrix[4][4] = shape[1][3]; matrix[4][5] = shape[1][2];
 matrix[4][7] = shape[2][3]; matrix[4][8] = shape[2][2]; matrix[4][10] = shape[3][3]; matrix[4][11] = shape[3][2];

  matrix[5][0] = shape[0][3]; matrix[5][2] = shape[0][1]; matrix[5][3] = shape[1][3]; matrix[5][5] = shape[1][1];
 matrix[5][6] = shape[2][3]; matrix[5][8] = shape[2][1]; matrix[5][9] = shape[3][3]; matrix[5][11] = shape[3][1];





}

//////////////////////////////////////////////////////////
float Determinants(float arr [4][4] , int m, int n)
{
float mat [3][3];
int k=0;

for(int i=0;i<4;i++)
if(i!=m)
{int l=0;
for(int j=0;j<4;j++)

if(j!=n)
{
mat[k][l] = arr[i][j];
l++;
}
k++;
}


float det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1]) - mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0]) +  mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);

return det;
}


///////////////////////////////////////////////////////////////////////////
void Matrix_Multiply(float  C[6][6], float B[6][12],float mat[6][12]){

for(int i=0;i<6;i++)

for(int j=0;j<12;j++)
{ float sum = 0;
for(int k=0;k<6;k++)
{
sum+=C[i][k]*B[k][j];
}
mat[i][j] = sum;
}

}

/////////////////////////////////////////////////////////////////////////////////
void Trans_Matrix_Multiply(float  B[6][12], float K[6][12],float mat[12][12]){



for(int i=0;i<12;i++)

for(int j=0;j<12;j++)
{ float sum = 0;
for(int k=0;k<6;k++)
{
sum+=B[k][i]*K[k][j];
}
mat[i][j] = sum;
}


}
