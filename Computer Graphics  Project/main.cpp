#include<GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include<math.h>
#include<bits/stdc++.h>
#include "include/BmpLoader.h"
#include "src/BmpLoader.cpp"
#include <GL/gl.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>


using namespace std;

const double PI = 3.14159265389;

double windowHeight=1920, windowWidth=1080;
unsigned int ID_[30] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
unsigned int sea_texture = ID_[9];
GLfloat angle = 0;
GLUquadricObj *qobj;
double Txval=0,Tyval=0,Tzval=0;

const int L=3;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 20;
GLfloat ctrlpoints[L+1][3] =
{
    {0.65,1.275,0.0},{2.45,1.25,0.0},{3.25,0.125,0.0}
};

long long nCr(int n, int r)
{
    if(r > n / 2)
        r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}


void bottleBezier(float col_a, float col_b, float col_c)
{
    GLfloat mat_ambient[] = { col_a*0.7, col_b*0.7, col_c*0.7, 1.0 };
    GLfloat mat_diffuse[] = {  col_a, col_b, col_c, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = {60};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);


    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}



static GLfloat v_pyramid[5][3] =
{
    {0.0, 0.0, 0.0},  //point index 0
    {0.0, 0.0, 2.0},  //point index 1
    {2.0, 0.0, 2.0},  //point index 2
    {2.0, 0.0, 0.0},  //point index 3
    {1.0, 4.0, 1.0}   //point index 4
};

static GLubyte p_Indices[4][3] =
{
    {4, 1, 2}, // indices for drawing the triangle plane 1
    {4, 2, 3}, // indices for drawing the triangle plane 2
    {4, 3, 0}, // indices for drawing the triangle plane 3
    {4, 0, 1}  // indices for drawing the triangle plane 4
};

GLfloat alpha = 0.0, theta = 0.0, axis_x=0.0, axis_y=0.0, fan_rot=360.0, f_x=0.0, f_y = 0.0, rot = 0.0, d_x=0.0, d_y=0.0;


static GLubyte quadIndices[1][4] =
{
    {0, 3, 2, 1}
};  // indeces for drawing the quad plane

static GLfloat v_cube[8][3] =
{
    /**
    {0.0, 0.0, 0.0},  //point index 0
    {2.0, 0.0, 0.0},  //point index 1
    {2.0, 0.0, 2.0},  //point index 2
    {0.0, 0.0, 2.0},  //point index 3
    {0.0, 2.0, 0.0},  //point index 4
    {2.0, 2.0, 0.0},  //point index 5
    {2.0, 2.0, 2.0},  //point index 6
    {0.0, 2.0, 2.0}   //point index 7
    **/

    {0,0,0},
    {0,0,2},
    {0,2,0},
    {0,2,2},

    {2,0,0},
    {2,0,2},
    {2,2,0},
    {2,2,2}
};

static GLubyte hexIndices[6][4] =
{
    /**
    {0, 1, 5, 4}, // indices for drawing the triangle plane 1
    {1, 2, 6, 5}, // indices for drawing the triangle plane 2
    {2, 3, 7, 6}, // indices for drawing the triangle plane 3
    {3, 0, 4, 7},  // indices for drawing the triangle plane 4
    {0, 1, 2, 3}, // indices for drawing the triangle plane 5
    {4, 5, 6, 7}  // indices for drawing the triangle plane 6
    **/

    /**
    {4,5,6,7},
    {1,2,6,5},
    {3,2,6,7},
    {0,3,7,4},
    {0,1,2,3},
    {0,1,5,4}**/

    {3,1,5,7},  //front
    {6,4,0,2},  //back
    {2,3,7,6},  //top
    {1,0,4,5},  //bottom
    {7,5,4,6},  //right
    {2,0,1,3}


};

GLboolean bRotate = false, uRotate = false, fanRotate = false, dRotate = false, ddRotate = false;
static GLfloat colors[6][3] =
{
    {0.0, 0.0, 1.0},  //color for point index 0
    {0.5, 0.0, 1.0},  //color for point index 1
    {0.0, 1.0, 0.0},  //color for point index 2
    {0.0, 1.0, 1.0},  //color for point index 3
    {0.8, 0.0, 0.0},  //color for point index 4
    {0.5, 1.0, 0.5}  //color for point index 5
};


static void getNormal_3p(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

static void getNormal3p(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void solid(float col_a,float col_b,float col_c)
{

    GLfloat mat_ambient[] = { col_a*0.7, col_b*0.7, col_c*0.7, 1.0 };
    GLfloat mat_diffuse[] = {  col_a, col_b, col_c, 1.0 };
    GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_shininess[] = {0};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glutSolidSphere (3.0, 20, 16);

}


void drawpyramid(float col_a,float col_b,float col_c)
{

    GLfloat mat_ambient[] = { col_a*0.7, col_b*0.7, col_c*0.7, 1.0 };
    GLfloat mat_diffuse[] = {  col_a, col_b, col_c, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = {60};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);


    ///glColor3f(col_a,col_b,col_c);
    glBegin(GL_TRIANGLES);

    for (GLint i = 0; i <4; i++)
    {
        //glColor3f(colors[i][0],colors[i][1],colors[i][2]);
        getNormal_3p(v_pyramid[p_Indices[i][0]][0], v_pyramid[p_Indices[i][0]][1], v_pyramid[p_Indices[i][0]][2],
                     v_pyramid[p_Indices[i][1]][0], v_pyramid[p_Indices[i][1]][1], v_pyramid[p_Indices[i][1]][2],
                     v_pyramid[p_Indices[i][2]][0], v_pyramid[p_Indices[i][2]][1], v_pyramid[p_Indices[i][2]][2]);

        //glColor3fv(&colors[i][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][0]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][1]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][2]][0]);
    }
    glEnd();

    glBegin(GL_QUADS);

    for (GLint i = 0; i <1; i++)
    {
        //glColor3f(colors[4][0],colors[4][1],colors[4][2]);
        getNormal_3p(v_pyramid[quadIndices[i][0]][0], v_pyramid[quadIndices[i][0]][1], v_pyramid[quadIndices[i][0]][2],
                     v_pyramid[quadIndices[i][1]][0], v_pyramid[quadIndices[i][1]][1], v_pyramid[quadIndices[i][1]][2],
                     v_pyramid[quadIndices[i][2]][0], v_pyramid[quadIndices[i][2]][1], v_pyramid[quadIndices[i][2]][2]);

        glVertex3fv(&v_pyramid[quadIndices[i][0]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][1]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][2]][0]);
        glVertex3fv(&v_pyramid[quadIndices[i][3]][0]);
    }
    glEnd();

}


GLboolean AL = false, DL = false, SL = false, L0=false, L1 = false, L2 = false,L3= false,L4=false,L5=false;



void drawCube(float col_a,float col_b,float col_c)
{

    GLfloat mat_ambient[] = { col_a*0.7, col_b*0.7, col_c*0.7, 1.0 };
    GLfloat mat_diffuse[] = {  col_a, col_b, col_c, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = {60};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    ///glColor3f(col_a,col_b,col_c);
    glBegin(GL_QUADS);

    for (GLint i = 0; i <6; i++)
    {
        //        glColor3f(colors[i][0],colors[i][1],colors[i][2]);
        getNormal3p(v_cube[hexIndices[i][0]][0], v_cube[hexIndices[i][0]][1], v_cube[hexIndices[i][0]][2],
                    v_cube[hexIndices[i][1]][0], v_cube[hexIndices[i][1]][1], v_cube[hexIndices[i][1]][2],
                    v_cube[hexIndices[i][2]][0], v_cube[hexIndices[i][2]][1], v_cube[hexIndices[i][2]][2]);

        glVertex3fv(&v_cube[hexIndices[i][0]][0]);
        glTexCoord2f(1,1);
        glVertex3fv(&v_cube[hexIndices[i][1]][0]);
        glTexCoord2f(1,0);
        glVertex3fv(&v_cube[hexIndices[i][2]][0]);
        glTexCoord2f(0,0);
        glVertex3fv(&v_cube[hexIndices[i][3]][0]);
        glTexCoord2f(0,1);
    }
    glEnd();
    //glutSolidSphere (3.0, 20, 16);

}
GLfloat  eye_x = 3, eye_y = 10, eye_z = 100, look_x=5, look_y=0.0;
/*void ownTranslatef(GLfloat dx, GLfloat dy, GLfloat dz)
{

    GLfloat m[16];

    m[0] = 1;
    m[4] = 0;
    m[8] = 0;
    m[12] = dx;
    m[1] = 0;
    m[5] = 1;
    m[9] = 0;
    m[13] = dy;
    m[2] = 0;
    m[6] = 0;
    m[10] = 1;
    m[14] = dz;
    m[3] = 0;
    m[7] = 0;
    m[11] = 0;
    m[15] = 1;

    glMatrixMode(GL_MODELVIEW);
    glMultMatrixf(m);
}
*/

double zx1 = 0.0,zx2 = 0.0,rotate1 = 0.0,rotate2 = 0.0,move_sun=0.0,move_moon = -3000.0;
void light()
{
    ///0th light(Spotlight)---headlight 01
    GLfloat no_light[] = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat light_ambient[]  = {1.0, 1.0, 1.0, 1.0};
    GLfloat light_diffuse[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = {31+27-3+0.5-2+0.5-2.5+6-1+18,13,10+4+20-4-29+5+173.5-1000+zx2-50,1};

    if(AL)
        glLightfv( GL_LIGHT0, GL_AMBIENT, light_ambient);
    else
        glLightfv( GL_LIGHT0, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    else
        glLightfv( GL_LIGHT0, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT0, GL_SPECULAR, light_specular);
    else
        glLightfv( GL_LIGHT0, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT0, GL_POSITION, light_position);


    GLfloat spot_direction[] = {0,0,1};
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spot_direction);
    glLightf( GL_LIGHT0, GL_SPOT_CUTOFF, 2);


    ///1st light(Spotlight)---headlight 02

    GLfloat light_position_1[] = {31+27-3+0.5-2+0.5-2.5+6-1,13,10+4+20-4-29+5+1000+zx1+50,1};

    if(AL)
        glLightfv( GL_LIGHT1, GL_AMBIENT, light_ambient);
    else
        glLightfv( GL_LIGHT1, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    else
        glLightfv( GL_LIGHT1, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT1, GL_SPECULAR, light_specular);
    else
        glLightfv( GL_LIGHT1, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT1, GL_POSITION, light_position_1);


    GLfloat spot_direction_1[] = {0,0,-1};
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, spot_direction_1);
    glLightf( GL_LIGHT1, GL_SPOT_CUTOFF, 2);

//    ///......................2nd light.........................GREEN
//    GLfloat light_ambient2[]  = {0.0, 1.0, 0.0, 1.0};
//    GLfloat light_diffuse2[]  = { 0.0, 1.0, 0.0, 1.0 };
//    GLfloat light_position2[] = { 10, 0, 80, 1.0 };
//
//    if(AL)
//        glLightfv( GL_LIGHT1, GL_AMBIENT, light_ambient2);
//    else
//        glLightfv( GL_LIGHT1, GL_AMBIENT, no_light);
//
//    if(DL)
//        glLightfv( GL_LIGHT1, GL_DIFFUSE, light_diffuse2);
//    else
//        glLightfv( GL_LIGHT1, GL_DIFFUSE, no_light);
//
//    if(SL)
//        glLightfv( GL_LIGHT1, GL_SPECULAR, light_specular);
//    else
//        glLightfv( GL_LIGHT1, GL_SPECULAR, no_light);
//
//    glLightfv( GL_LIGHT1, GL_POSITION, light_position2);


///......................2nd light.........................///
    GLfloat light_ambient2[]  = {1, 1, 1, 1.0};
    GLfloat light_diffuse2[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position2[] = {-1000,130,1000+move_sun, 0.0 };

    if(AL)
        glLightfv( GL_LIGHT2, GL_AMBIENT, light_ambient2);
    else
        glLightfv( GL_LIGHT2, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
    else
        glLightfv( GL_LIGHT2, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT2, GL_SPECULAR, light_specular);
    else
        glLightfv( GL_LIGHT2, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT2, GL_POSITION, light_position2);


///......................3th light.........................///
    GLfloat light_ambient3[]  = {1.0, 1.0, 0.0, 1.0};
    GLfloat light_diffuse3[]  = { 1.0, 1.0, 0.0, 1.0 };
    GLfloat light_position3[] = { -1000,130,1000+move_sun, 0.0 };

    if(AL)
        glLightfv( GL_LIGHT3, GL_AMBIENT, light_ambient3);
    else
        glLightfv( GL_LIGHT3, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT3, GL_DIFFUSE, light_diffuse3);
    else
        glLightfv( GL_LIGHT3, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT3, GL_SPECULAR, light_specular);
    else
        glLightfv( GL_LIGHT3, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT3, GL_POSITION, light_position3);


    ///......................4th light.........................///
    GLfloat light_ambient4[]  = {0.7, 0.7, 0.7, 1.0};
    GLfloat light_diffuse4[]  = { 0.7, 0.7, 0.7, 1.0 };
    GLfloat light_position4[] = { -1000,130,1000+move_sun, 0.0 };

    if(AL)
        glLightfv( GL_LIGHT4, GL_AMBIENT, light_ambient4);
    else
        glLightfv( GL_LIGHT4, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT4, GL_DIFFUSE, light_diffuse4);
    else
        glLightfv( GL_LIGHT4, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT4, GL_SPECULAR, light_specular);
    else
        glLightfv( GL_LIGHT4, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT4, GL_POSITION, light_position4);


    ///......................5th light.........................///
    GLfloat light_ambient5[]  = {0.1, 0.1, 0.1, 1.0};
    GLfloat light_diffuse5[]  = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat light_specular5[]  = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat light_position5[] = { 500,130,-1000+move_moon, 0.0 };

    if(AL)
        glLightfv( GL_LIGHT5, GL_AMBIENT, light_ambient5);
    else
        glLightfv( GL_LIGHT5, GL_AMBIENT, no_light);

    if(DL)
        glLightfv( GL_LIGHT5, GL_DIFFUSE, light_diffuse5);
    else
        glLightfv( GL_LIGHT5, GL_DIFFUSE, no_light);

    if(SL)
        glLightfv( GL_LIGHT5, GL_SPECULAR, light_specular5);
    else
        glLightfv( GL_LIGHT5, GL_SPECULAR, no_light);

    glLightfv( GL_LIGHT5, GL_POSITION, light_position5);

}



void building(void)
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();
    glTranslatef(-25+2,-3,-25);
    glScalef(0.05,10,50);
    drawCube(1,1,1);
    glPopMatrix();


}

void building2(void)
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();
    glTranslatef(-25+2,-3,-25+0.05);
    glScalef(37,10,0.05);
    drawCube(1,1,1);
    glPopMatrix();


}

void building3(void)
{

    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();
    glTranslatef(-25+2+98-0.05,-3,-25);
    glScalef(0.05,10,50);
    drawCube(1,1,1);
    glPopMatrix();

}

void railway_tunnel(void)
{
    /**glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(-25+2+74,-3,-25+0.05);
    glScalef(20,5.05,0.05);
    drawCube(0,0,0);
    glPopMatrix();
    **/

    glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(-25+2+74-4,-3,-400);
    glScalef(1,8.5,200);
    drawCube(1,1,1);
    glPopMatrix();

    glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(-25+2+74+30,-3,-400);
    glScalef(1,8.5,200);
    drawCube(1,1,1);
    glPopMatrix();

    glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(-25+2+74+20-7,-3,-400);
    glScalef(1,8.5,200);
    drawCube(1,1,1);
    glPopMatrix();

    glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(-25+2+74-4,-3+17,-400);
    glScalef(18,2,200);
    drawCube(1,1,1);
    glPopMatrix();

}

void sky()
{
    glViewport(0, 0, windowHeight, windowWidth);



    glPushMatrix();
    glTranslatef(-25-1000,-3+200,-1000);
    glScalef(2000,0.05,1000);
    drawCube(1,1,1);
    glPopMatrix();



    glPushMatrix();
    glTranslatef(-25-1000,-3-55,-1000);
    glScalef(0.05,130,1000);
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-25+1000,-3-55,-1000);
    glScalef(0.05,130,1000);
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-25-1000,-3-55,-1000);
    glScalef(2000,130,0.05);
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-25-1000,-3-55,1000);
    glScalef(2000,130,0.05);
    drawCube(1,1,1);
    glPopMatrix();







}



void rail (int xx, int zz, int xrotate)
{





    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27-2.5,-1,10+zz);    ///1
    glRotatef(xrotate,1,0,0);
    glRotatef(180,0,1,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();


    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27-2.5,-1,10+4+zz);    ///2
    glRotatef(xrotate,1,0,0);
    glRotatef(180,0,1,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();

    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27.5,-1,10+zz);         ///3
    glRotatef(xrotate,1,0,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();


    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27.5,-1,10+4+zz);          ///4
    glRotatef(xrotate,1,0,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();



    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27-2.5,-1,10+20+zz);
    glRotatef(xrotate,1,0,0);
    glRotatef(180,0,1,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();


    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27-2.5,-1,10+4+20+zz);
    glRotatef(xrotate,1,0,0);
    glRotatef(180,0,1,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();

    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27.5,-1,10+20+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();


    glPushMatrix();
    //glColor3f(0,0,0);
    glTranslatef(xx+31+27.5,-1,10+4+20+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(0.5,1.5,1.5);
    bottleBezier(0,0,0);
    glPopMatrix();



    /**** wheel axis ***/

    glPushMatrix();
    glTranslatef(xx+31+27-3,-1,10+4+20-4-20+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(2.2,0.2,0.2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3,-1,10+4+20-20+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(2.2,0.2,0.2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3,-1,10+4+20-4+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(2.2,0.2,0.2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3,-1,10+4+20+zz);
    glRotatef(xrotate,1,0,0);
    glScalef(2.2,0.2,0.2);
    drawCube(0,0,0);
    glPopMatrix();

    /*** support ***/
    ///1
    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5,-1,10+4+20-4-20+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+3,-1,10+4+20-4-20+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-1,-1+2,10+4+20-4-20+zz);
    glScalef(2.5,0.5,0.5);
    drawCube(0,0,0);
    glPopMatrix();


    ///2
    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5,-1,10+4+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+3,-1,10+4+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-1,-1+2,10+4+zz);
    glScalef(2.5,0.5,0.5);
    drawCube(0,0,0);
    glPopMatrix();


    ///3
    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5,-1,10+4+20-4+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+3,-1,10+4+20-4+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-1,-1+2,10+4+20-4+zz);
    glScalef(2.5,0.5,0.5);
    drawCube(0,0,0);
    glPopMatrix();


    ///4
    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5,-1,10+4+20+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+3,-1,10+4+20+zz);
    glScalef(0.3,1.1,0.3);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-1,-1+2,10+4+20+zz);
    glScalef(2.5,0.5,0.5);
    drawCube(0,0,0);
    glPopMatrix();


    /*** body ***/

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5,2,10+4+20-4-20-4+zz);
    glScalef(5,6,16);
    drawCube(0,0,1);
    glPopMatrix();

    /** door1 **/

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05,2,10+4+20-4-20-2+zz);
    glScalef(0.05,5,2);
    drawCube(1,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2+10,10+4+20-4-20-2+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+4-0.49+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();


    /** door2 **/
    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,5,2);
    drawCube(1,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2+10,10+4+20-4-20-2+24+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1,2,10+4+20-4-20-2+4-0.49+24+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();


    /** door3 **/

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10,2,10+4+20-4-20-2+zz);
    glScalef(0.05,5,2);
    drawCube(1,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2+10,10+4+20-4-20-2+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+4-0.49+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();


    /** door4 **/

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,5,2);
    drawCube(1,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+24+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2+10,10+4+20-4-20-2+24+zz);
    glScalef(0.05,0.25,2);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.1+10+0.06,2,10+4+20-4-20-2+4-0.49+24+zz);
    glScalef(0.05,5,0.25);
    drawCube(0,0,0);
    glPopMatrix();

    /** window **/
    for(int i=9; i<=33; i+=4)
    {
        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10,7,i+zz);
        glScalef(0.05,1,1);
        drawCube(1,1,0);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10+0.05,7,i+zz);
        glScalef(0.05,1,0.25);
        drawCube(0,0,0);
        glPopMatrix();
    }
    for(int i=13; i<=29; i+=4)
    {
        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10+0.05,7,i+2-0.25+zz);
        glScalef(0.05,1,0.25);
        drawCube(0,0,0);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05+10+0.05,7+2,i+zz);
        glScalef(0.05,0.25,1);
        drawCube(0,0,0);
        glPopMatrix();
    }


    for(int i=9; i<=33; i+=4)
    {
        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.05,7,i+2+zz);
        glRotatef(180,0,1,0);
        glScalef(0.05,1,1);
        drawCube(1,1,0);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.12,7,i+2.1+zz);
        glRotatef(180,0,1,0);
        glScalef(0.05,1,0.25);
        drawCube(0,0,0);
        glPopMatrix();
    }
    for(int i=13; i<=29; i+=4)
    {
        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.12,7,i+2-0.25+zz);
        glRotatef(180,0,1,0);
        glScalef(0.05,1,0.25);
        drawCube(0,0,0);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(xx+31+27-3+0.5-2+0.5-2.5-0.12,7+2,i+2.1+zz);
        glRotatef(180,0,1,0);
        glScalef(0.05,0.25,1);
        drawCube(0,0,0);
        glPopMatrix();
    }

    /*** connection between two trains ***/

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5+5,2+4,10+4+20-4-20-4-1.5+zz);
    glScalef(0.50,0.50,0.75);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5+5,2+3.5,10+4+20-4-20-4-1.5+zz);
    glScalef(1,1,0.25);
    drawCube(0,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5+5+1,2+4,10+4+20-4-20-4-1.5+35+zz);
    glRotatef(180,0,1,0);
    glScalef(0.50,0.50,0.75);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(xx+31+27-3+0.5-2+0.5-2.5+5+1,2+3.5,10+4+20-4-20-4-1.5+35+zz);
    glRotatef(180,0,1,0);
    glScalef(1,1,0.25);
    drawCube(0,0,0);
    glPopMatrix();






}

void sea(void)
{
    glPushMatrix();
    glPushMatrix();
    glTranslatef(-1000,-50,-1000);
    glScalef(500,0.05,1000);
    drawCube(1,1,1);
    glPopMatrix();
    glPopMatrix();
}


void railline()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();
    glTranslatef(27+28,-3,-1000);
    glScalef(0.1,0.08,1000);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(31+28,-3,-1000);
    glScalef(0.1,0.08,1000);
    drawCube(0,0,0);
    glPopMatrix();

    for(double z=-20; z<=1500; z+=1.5)
    {
        glPushMatrix();
        glTranslatef(29-2.5+28,-3,-980+z);
        glScalef(2.5,0.08,0.4);
        drawCube(0.5,0.2,0.0);
        glPopMatrix();
    }


    int xx = 18;

    glPushMatrix();
    glTranslatef(27+28+xx,-3,-1000);
    glScalef(0.1,0.08,1000);
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(31+28+xx,-3,-1000);
    glScalef(0.1,0.08,1000);
    drawCube(0,0,0);
    glPopMatrix();

    for(double z=-20; z<=1500; z+=1.5)
    {
        glPushMatrix();
        glTranslatef(29-2.5+28+xx,-3,-980+z);
        glScalef(2.5,0.08,0.4);
        drawCube(0.5,0.2,0.0);
        glPopMatrix();
    }

}


void toy(void)
{
    int n=4;
    GLfloat rotation = 0;
    for(int i=0; i<n; i++)
    {

        glPushMatrix();
        glRotatef(rotation,0,1,0);
        glTranslatef(0,2,0);
        glScalef(2,2,2);
        drawCube(0,0,0);
        glPopMatrix();
        rotation+= (GLfloat)360.0/n;
    }

}
void road(void)
{

    glViewport(0, 0, windowHeight, windowWidth);



    glEnable(GL_TEXTURE_2D);
    glPushMatrix();

    glBindTexture(GL_TEXTURE_2D,ID_[7]);
    glPushMatrix();
    glTranslatef(-25,-3,-10-16-1000);
    glScalef(5.5,0.07,1000);
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-25+72-8,-3,-10-16-1000);
    glScalef(4,0.07,1000);
    drawCube(1,1,1);
    glPopMatrix();

    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();

    glBindTexture(GL_TEXTURE_2D,ID_[8]);

    glPushMatrix();
    glTranslatef(-25+11,-3,-10-8-7);
    glScalef(26.5,0.07,3.8);
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-25+11,-3,-10-8-7+62);
    glScalef(26.5,0.07,3.8);
    drawCube(1,1,1);
    glPopMatrix();

    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


}
void deyal()
{
    glViewport(0, 0, windowHeight, windowWidth);

    ///glPushMatrix(); //...........begin of walls..............

    //ground
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID_[1]);
    glPushMatrix();
    glTranslatef(-100+50,-3,-1000);
    glScalef(1000,0.05,1000);
    drawCube(1,1,1);
    glPopMatrix();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID_[0]);
    //wall1
    glPushMatrix();
    glTranslatef(-14,-3,-10-8);
    glScalef(26,2,0.5);
    drawCube(1,1,1);
    glPopMatrix();

    //wall2_1
    glPushMatrix();
    glTranslatef(-14,-3,-10-8);
    glScalef(0.5,2,13);
    drawCube(1,1,1);
    glPopMatrix();

    //wall2_2
    glPushMatrix();
    glTranslatef(-14,-3,12);
    glScalef(0.5,2,13);
    drawCube(1,1,1);
    glPopMatrix();

    //wall3_1
    glPushMatrix();
    glTranslatef(-14+60-8,-3,-10-8);
    glScalef(0.5,2,13);
    drawCube(1,1,1);
    glPopMatrix();

    //wall3_2
    glPushMatrix();
    glTranslatef(-14+60-8,-3,12);
    glScalef(0.5,2,13);
    drawCube(1,1,1);
    glPopMatrix();

    //wall4
    glPushMatrix();
    glTranslatef(-14,-3,38);
    glScalef(26,2,0.5);
    drawCube(1,1,1);
    glPopMatrix();



    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void tree()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix(); //...........begin of tree..............
    glTranslatef(0, 0.2, 0);

    //up

    glPushMatrix();
    glColor4f(0,1,0,100);
    glTranslatef(0.2,1.5,0);
    //drawpyramid(0.1,0.5,0.0);
    glutSolidSphere(2,60,60);
    glPopMatrix();


    //down
    glPushMatrix();
    glTranslatef(0.5,-3,0.5);
    glScalef(0.5,2.5,0.5);
    drawCube(0.5,0.2,0.0);
    glPopMatrix();

    glPopMatrix();
}

void array_of_tree()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();

    //back
    glPushMatrix();
    glTranslatef(0,0,-9);
    glScalef(0.8,1,0.8);
    for(int i = -8; i<28; i+=4)
    {
        glPushMatrix();
        glTranslatef(i,0,0);
        tree();
        glPopMatrix();
    }
    glPopMatrix();


    //front
    glPushMatrix();
    glTranslatef(0,0,18);
    glScalef(0.8,1,0.8);
    for(int i = -8; i<28; i+=4)
    {
        glPushMatrix();
        glTranslatef(i,0,0);
        tree();
        glPopMatrix();
    }
    glPopMatrix();


    //left
    glPushMatrix();
    glTranslatef(-7,0,0);
    glScalef(0.8,1,0.8);
    for(int i = -8; i<20; i+=5)
    {
        glPushMatrix();
        glTranslatef(0,0,i);
        tree();
        glPopMatrix();
    }
    glPopMatrix();


    //right
    glPushMatrix();
    glTranslatef(20,0,0);
    glScalef(0.8,1,0.8);
    for(int i = -8; i<20; i+=5)
    {
        glPushMatrix();
        glTranslatef(0,0,i);
        tree();
        glPopMatrix();
    }
    glPopMatrix();


    glPopMatrix();
}

void ladder()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();

    glPushMatrix();
    glTranslatef(0,-3,0);
    glScalef(0.2,2.5,0.2);
    drawCube(0.2,0.1,0.4);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2,-3,0);
    glScalef(0.2,2.5,0.2);
    drawCube(0.2,0.1,0.4);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(0,-3,0);
    glScalef(1.2,0.1,0.2);
    for(int i = 3; i<44; i+=10)
    {
        glPushMatrix();
        glTranslatef(0,i,0);
        drawCube(0.4,0.1,0.1);
        glPopMatrix();
    }
    glPopMatrix();

    glPopMatrix();

}

void call_ladder()
{

    glPushMatrix();
    ladder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,-6);
    ladder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,2,-2);
    glScalef(1,1,1.2);
    glRotated(90,1,0,0);
    ladder();
    glPopMatrix();
}

void sliding_ladder()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();

    glPushMatrix();
    ladder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,1.8,0);
    glScalef(1.2,0.1,1);
    drawCube(0.4,0.2,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,2,2);
    glRotatef(140,1,0,0);
    glScalef(1.2,3,0.1);
    drawCube(0.4,0.2,0.1);
    glPopMatrix();

    glPopMatrix();
}

void chorki_stand(int z)
{
    ///front triangle
    glPushMatrix();
    glPushMatrix();
    glTranslatef(20,-3,0+z);
    glRotatef(15, 0, 0, 1);
    glScalef(0.8,10,0.8);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1,-3,0+z);
    glRotatef(-15, 0, 0, 1);
    glScalef(0.8,10,0.8);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPopMatrix();
}

void dolna_stand_func()
{
    ///front triangle
    glPushMatrix();

    glPushMatrix();
    glRotatef(15, 0, 0, 1);
    glTranslatef(0,-3,0);
    glScalef(0.2,1.3,0.2);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-15, 0, 0, 1);
    glTranslatef(-1,-3,0);
    glScalef(0.2,1.3,0.2);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPopMatrix();

    ///back triangle
    glPushMatrix();
    glTranslatef(0,0,-5);

    glPushMatrix();
    glRotatef(15, 0, 0, 1);
    glTranslatef(0,-3,0);
    glScalef(0.2,1.3,0.2);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-15, 0, 0, 1);
    glTranslatef(-1,-3,0);
    glScalef(0.2,1.3,0.2);
    drawpyramid(1.0, 0, 0);
    glPopMatrix();

    glPopMatrix();

    ///upper rod
    glPushMatrix();
    glTranslatef(-0.4,1.9,-4.8);
    glScalef(0.1, 0.1, 2.5);
    drawCube(1.0, 0, 0);
    glPopMatrix();

}

void dolna_sitting_place_func()
{
    //connector
    glPushMatrix();
    glTranslatef(-1.5,0,-0.4);
    glScalef(1,0.5,0.5);
    drawpyramid(0.5, 0, 0);
    glPopMatrix();

    //sitting place
    glPushMatrix();
    glRotatef( rot, d_x, d_y, 0.0 );

    //hanging stand
    glPushMatrix();

    glPushMatrix();
    glRotatef(10, 0, 0, 1);
    glTranslatef(0,-3,0);
    glScalef(0.1,0.8,0.1);
    drawpyramid(0.5, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-10, 0, 0, 1);
    glTranslatef(-1,-3,0);
    glScalef(0.1,0.8,0.1);
    drawpyramid(0.5, 0, 0);
    glPopMatrix();

    glPopMatrix();

    //flat
    glPushMatrix();
    glTranslatef(-2,-3,-1);
    glScalef(1.5,0.2,1);
    drawCube(0.4,0,0);
    glPopMatrix();

    glPopMatrix();
}

void call_dolna_func()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glTranslated(2,0,0);

    glPushMatrix();//...............dolna............

    glPushMatrix();
    glTranslatef(1,0,-0.3);
    glRotatef(90,0,1,0);
    dolna_stand_func();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,0);
    glScalef(0.5,1,1);
    dolna_sitting_place_func();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-2,0,0);
    glScalef(0.5,1,1);
    dolna_sitting_place_func();
    glPopMatrix();

    glPopMatrix();//.................end dolna............

}









void nagordola_stand()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();    //............begin of celling fan................

    //middle rod
    glPushMatrix();
    glTranslatef(-0.3, -2, -0.3);
    glScalef(0.3,2,0.3);
    drawCube(0.3,0.5,0.8);
    glPopMatrix();

//middle bar
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-1, -0.2, -1);
    glScalef(1,0.2,1);
    drawCube(0.5, 0.2, 1.0);
    glPopMatrix();


    //left-right wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-5.5, 0, -0.5);
    glScalef(6,0.05,0.5);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    //front-back wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-0.5, 0, -5.5);
    glScalef(0.5,0.05,6);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    glPopMatrix(); //............end of celling fan................

}


void nagordola_stand2()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();    //............begin of celling fan................

    //middle rod
    glPushMatrix();
    glTranslatef(-0.3, -2, -0.3);
    glScalef(0.3,2,0.3);
    drawCube(0.3,0.5,0.8);
    glPopMatrix();

    //wing connectors
    glPushMatrix(); //up
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-0.3, -4, -5.3);
    glScalef(0.3,2,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix(); //down
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-0.3, -4, 5.5);
    glScalef(0.3,2,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix(); //left
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-5.3, -4, -0.3);
    glScalef(0.3,2,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix(); //right
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(5.7, -4, -0.3);
    glScalef(0.3,2,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

//middle bar
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-1, -0.2, -1);
    glScalef(1,0.2,1);
    drawCube(0.5, 0.2, 1.0);
    glPopMatrix();


    //left-right wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-5.5, 0, -0.5);
    glScalef(6,0.05,0.5);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    //front-back wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-0.5, 0, -5.5);
    glScalef(0.5,0.05,6);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    glPopMatrix(); //............end of celling fan................

}


void nagordola_func()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();    //............begin of nagordola................

    //ones side
    glPushMatrix();
    glTranslatef(0,5,0);
    glRotatef(90,1,0,0);
    nagordola_stand();
    glPopMatrix();

    //other side
    glPushMatrix();
    glTranslatef(0,5,4);
    glRotatef(90,1,0,0);
    nagordola_stand2();
    glPopMatrix();

    //stand1
    glPushMatrix();
    glTranslatef(-1,-2.9,-2);
    glScalef(1,2.2,0.5);
    drawpyramid(0.2,0.2,0.5);
    glPopMatrix();

    //stand2
    glPushMatrix();
    glTranslatef(-1,-2.9,5);
    glScalef(1,2.2,0.5);
    drawpyramid(0.2,0.2,0.5);
    glPopMatrix();

//    //bike1
//    glPushMatrix();
//    glRotatef( fan_rot, f_x, f_y, 0.0 );
//     //glTranslatef(-1,-2.9,5);
//    bike_help();
//    glPopMatrix();


    glPopMatrix(); //............end of nagordola......
}







void dheki_sitting_place()
{
    //one side
    glPushMatrix();
    //handle
    glPushMatrix();
    glTranslatef(-1.2,0,0);

    glPushMatrix();
    glTranslatef(2.5,-2,0);
    glScalef(0.2, 0.7, 0.2);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-2,1);
    glScalef(0.2, 0.7, 0.2);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-1,0);
    glScalef(0.2, 0.2, 0.7);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPopMatrix();


    //stand
    glPushMatrix();
    glTranslatef(-2.5,-3,1.4);
    glScalef(1, 0.7, 1);
    glRotatef(180,1,0,0);

    glPushMatrix();
    glTranslatef(2.5,-2,0);
    glScalef(0.2, 1.2, 0.2);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-2,1);
    glScalef(0.2, 1.2, 0.2);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-0.5,0);
    glScalef(0.2, 0.2, 0.7);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPopMatrix();

    //flat
    glPushMatrix();
    glTranslatef(-0.1,-2.1,-0.1);
    glScalef(1, 0.3, 0.8);
    drawCube(0.6,0.7, 0.1);
    glPopMatrix();

    glPopMatrix();
}

void dheki_side()
{
    glPushMatrix();
    glTranslatef(2.5,-3,0);
    glScalef(0.2, 1.5, 0.3);
    drawCube(0.5,0.5, 0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-3,1);
    glScalef(0.2, 1.5, 0.3);
    drawCube(0.5,0.5, 0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2.5,-1,0);
    glScalef(0.2, 0.3, 0.6);
    drawCube(0.5,0.5, 0.5);
    glPopMatrix();

}

void call_dheki_func()
{
    glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();//..........start.......

    //middle bar
    glPushMatrix();
    glTranslatef(-0.7,1.8,2.5);
    glRotatef(90,0,1,0);

    glPushMatrix();
    glTranslatef(2.5,-2.2,0.2);
    glScalef(0.2, 0.2, 1.2);
    drawCube(0.5,0.2, 0.2);
    glPopMatrix();
    glPopMatrix();

    //side stand
    glPushMatrix();
    glTranslatef(-1, 0, -1);
    dheki_side();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-3, 0, -1);
    dheki_side();
    glPopMatrix();


    glPushMatrix();//..........begin moving part........
    glRotatef( rot, d_x, d_y, 0.0 );

    glPushMatrix();
    glTranslatef(0,2,3.2);
    glRotatef(90,0,1,0);

    glPushMatrix();
    glTranslatef(0,-2,0);
    glScalef(3, 0.2, 0.2);
    drawCube(0.5,0.2, 0.0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,-2,1);
    glScalef(3, 0.2, 0.2);
    drawCube(0.5,0.2, 0.0);
    glPopMatrix();

    glPushMatrix();
    dheki_sitting_place();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.5,0,1.4);
    glRotatef(180,0,1,0);
    dheki_sitting_place();
    glPopMatrix();

    glPopMatrix();

    glPopMatrix();

    glPopMatrix();//............end.........
}


void bike_help()
{
    //head
    glPushMatrix();
    glTranslatef(0,0.1,0);
    glScalef(1,0.8,1);
    drawpyramid(0.5,0.1,0.1);
    glPopMatrix();

    //sitting place
    glPushMatrix();
    glScalef(2,0.5,1);
    drawCube(0.5,0.2,0.5);
    glPopMatrix();

    //handle
    glPushMatrix();
    glTranslatef(0.8, 1.7, 0);
    glScalef(0.2,0.2,1);
    drawCube(0.5,0.2,0.5);
    glPopMatrix();

    //feet-holder
    glPushMatrix();
    glTranslatef(0.8, 0.1, -0.5);
    glScalef(0.2,0.2,1.5);
    drawCube(0.5,0.0,0.1);
    glPopMatrix();

    //back-wall
    glPushMatrix();
    glTranslatef(4, 0, 0);
    glScalef(0.2,1.5,1);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

    //roof
    glPushMatrix();
    glTranslatef(0, 2.9, 0);
    glScalef(2,0.2,1);
    drawpyramid(0.5,0.0,0.0);
    glPopMatrix();
}

void bike_func()
{
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();    //............begin of bike................

    //roof
    glPushMatrix();
    glTranslatef(-8, 1, -7);
    glScalef(8,0.5,8);
    drawpyramid(0.2,0.2,0.5);
    glPopMatrix();

    //ceil
    glPushMatrix();
    glTranslatef(-9, -9.2, -9);
    glScalef(9,0.2,9);
    drawpyramid(0.2,0.5,0.5);
    glPopMatrix();

    //string1
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-5, -6, 0);
    glScalef(0.2,3,0.2);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

    //string2
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(5.5, -6, 0);
    glScalef(0.2,3,0.2);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

    //string3
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(0, -6, -5);
    glScalef(0.2,3,0.2);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

    //string4
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(0, -6, 5.5);
    glScalef(0.2,3,0.2);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

    //sitting place front
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-1.95,-9,4.7);
    bike_help();
    glPopMatrix();


    //sitting place back
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(2,-9,-3.8);
    glRotatef(180,0,1,0);
    bike_help();
    glPopMatrix();

    //sitting place right
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(4.8,-9,2);
    glRotatef(90,0,1,0);
    bike_help();
    glPopMatrix();


    //sitting place left
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-3.8,-9,-2);
    glRotatef(270,0,1,0);
    bike_help();
    glPopMatrix();


    //middle rod
    glPushMatrix();
    glTranslatef(-0.3, -9.5, -0.3);
    glScalef(0.3,6,0.3);
    drawCube(0.2,0.2,0.5);
    glPopMatrix();

//middle bar
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-1, -0.2, -1);
    glScalef(1,0.2,1);
    drawCube(0.5, 0.2, 1.0);
    glPopMatrix();


    //left-right wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-5.5, 0, -0.5);
    glScalef(6,0.05,0.5);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    //front-back wing
    glPushMatrix();
    glRotatef( fan_rot, f_x, f_y, 0.0 );
    glTranslatef(-0.5, 0, -5.5);
    glScalef(0.5,0.05,6);
    drawCube(0.5, 0.5, 0.8);
    glPopMatrix();

    glPopMatrix(); //............end of celling fan................


}


///......................................................................

void window()
{

    glPushMatrix();
    glTranslatef(-4.9,2.5,-1.5);
    glScalef(0.35,2,0.1);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.9,2.5,-4.2);
    glScalef(0.35,2,0.1);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.9,6.1,-4.2);
    glScalef(0.35,0.2,1.45);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.9,2.3,-4.2);
    glScalef(0.35,0.2,1.45);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.4,5.4,-4.2);
    glScalef(0.05,0.05,1.45);
    drawCube(.2,.2,.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.4,4.4,-4.2);
    glScalef(0.05,0.05,1.45);
    drawCube(.2,.2,.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.4,3.4,-4.2);
    glScalef(0.05,0.05,1.45);
    drawCube(.2,.2,.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4.4,2.6,-2.9);
    glScalef(0.05,1.9,0.05);
    drawCube(.2,.2,.2);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-4.9,2.5,-4);
    glScalef(0.25,2,1.3);
    drawCube(1,1,1);
    glPopMatrix();

}

void window2()
{

    glPushMatrix();
    glTranslatef(5,2,-4.9);
    glScalef(0.1,2,0.35);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7.5,2,-4.9);
    glScalef(0.1,2,0.35);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,6,-4.9);
    glScalef(1.35,0.1,0.35);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,2,-4.9);
    glScalef(1.35,0.1,0.35);
    drawCube(0.52,0.37,0.26);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,3,-4.9);
    glScalef(1.35,0.05,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,4,-4.9);
    glScalef(1.35,0.05,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,5,-4.9);
    glScalef(1.35,0.05,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.25,2,-4.9);
    glScalef(0.07,2,0.3);
    drawCube(0.2,0.2,0.2);
    glPopMatrix();




    glPushMatrix();
    glTranslatef(5,2,-4.9);
    glScalef(1.3,2,0.25);
    drawCube(1,1,1);
    glPopMatrix();
}

void fanparts()
{
    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    ///glRectf(-1, -1, 1, 1);


    glTranslatef(2,6.5,0);
    glScalef(0.25,0.1,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,0);
    glScalef(1,0.001,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,0);
    glScalef(0.25,0.001,1);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(1,6.55,0);
    glScalef(1,0.001,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,-1);
    glScalef(0.25,0.001,1);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();
}
void building_help()
{
    GLfloat reduce = -0.5;
    glPushMatrix();
    glTranslatef(2+reduce,0,-5);
    glScalef(2,0.3,4);          ///khat
    drawCube(0.5,0.0,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2+reduce,0.6,-3);
    glScalef(0.9,0.1,1);        /// balish 1
    drawCube(1,0.8,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(4.1+reduce,0.6,-3);
    glScalef(0.9,0.1,1);        /// balish 2
    drawCube(1,0.8,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2+reduce,0.6,0);
    glScalef(1.9,0.1,1.2);        /// bedsheet
    drawCube(2/5,1/5,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2+reduce,0.6,-4);
    glScalef(2,0.6,0.1);        /// bed head
    drawCube(0.5,0.0,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0.6,0);
    glScalef(0.6,0.0,1.3);        /// carpet
    drawCube(1,0.0,1);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(6.5,0.0,-4);
    glScalef(0.8,0.8,0.8);        /// drawer 1
    drawCube(0.5,0.0,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1,0.0,-4);
    glScalef(0.8,1.8,1);        /// drawer 2
    drawCube(0.5,0.0,0.1);
    glPopMatrix();


    GLfloat ex = 1.75;
    glPushMatrix();
    glTranslatef(-4,0.0,3-ex);
    glScalef(0.8,1.2,1.8);        /// drawer 3
    drawCube(0.5,0.0,0.1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4,12/5,3-ex);
    glScalef(0,1.5,1.5);        /// mirror
    drawCube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4,12/5 + 3,3-ex);  ///mirror
    glScalef(0.1,0.1,1.5);
    drawCube(1,1,153/255);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4,12/5 + 3,3-ex);   ///mirror
    glRotatef(90,1,0,0);
    glScalef(0.1,0.1,1.5);
    drawCube(1,1,153/255);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-4,12/5 + 3,6-0.2-ex);
    glRotatef(90,1,0,0);                     ///mirror
    glScalef(0.1,0.1,1.5);
    drawCube(1,1,153/255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-2,5,-3);
    glScalef(1.0,1,0.1);        /// picture
    drawCube(0.9,0.9,0.6);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1.25,6,-3);
    glScalef(0.2,0.2,0.12);        /// pic_head
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1.25,6-0.25,-3);
    glScalef(0.1,0.18,0.12);        /// pic_throat
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-1.75,6-0.95,-3);
    glScalef(0.7,0.45,0.12);        /// pic_body
    drawCube(0,0.0,0.0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(1,5+0.5,-3);
    glScalef(0.5,0.5,0.1);        /// clock
    drawCube(0.95,0.95,0.95);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(1+0.5-0.1,6-0.25+0.2,-3+0.21);
    glScalef(0.08,0.08,0);        /// clock_dot
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(1+0.5-0.1,6-0.25+0.2,-3+0.21);
    glScalef(0.25,0.05,0);        /// clock_hand1
    drawCube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(1+0.5-0.1,6-0.25+0.2,-3+0.21);
    glScalef(0.05,0.25,0);        /// clock_hand2
    drawCube(0,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-5,8,-5);     /// ceiling
    glScalef(7,0,7);
    drawCube(0.8,0.8,0.8);
    glPopMatrix();

    window();
    window2();


    ///number 1 cube

    ///The floor
    for(int i=-10; i<18; i++)
    {
        for(int j=-10; j<18; j++)
        {
            glPushMatrix();
            glTranslatef(0.5*j,0.5,0.5*i);
            glScalef(0.25,0.02,0.25);
            if((i+j)%2 == 0)
            {
                drawCube(1,1,1);
            }
            else
            {
                drawCube(0,0,0);
            }
            glPopMatrix();

        }

    }

    glPushMatrix();
    glTranslatef(-6,-1+8+0.05,-5-1);
    glScalef(8,0.5,8);
    drawCube(0.1,0.1,0.2);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(2,6.7,0);
    glScalef(0.05,0.3,0.05);
    drawCube(0.1,0.1,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    ///glRectf(-1, -1, 1, 1);


    glTranslatef(2,6.5,0);
    glScalef(0.25,0.1,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    ///glRotatef(angle,0,1,0);

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,0);
    glScalef(1,0.001,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,0);
    glScalef(0.25,0.001,1);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(1,6.55,0);
    glScalef(1,0.001,0.25);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    ///glRotatef(theta,0,1,0);
    glTranslatef(2,6.55,-1);
    glScalef(0.25,0.001,1);
    drawCube(0.1,0.2,0.2);
    glPopMatrix();
}

void building_wall()
{
    glPushMatrix();

    glPushMatrix();
    glTranslatef(-10*0.5,0.0,-10*0.5);
    glScalef(0.25,4,7);
    drawCube(0.8,0.8,0.8);                         ///wall 04
    glPopMatrix();

    //glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(18*0.5,0.0,-10*0.5);
    glScalef(0.25,4,7);
    drawCube(0.8,0.8,0.8);                       ///wall 03
    glPopMatrix();

    //glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(18*0.5,0.0,-10*0.5);
    glRotatef(-90,0,1,0);               ///wall 02
    glScalef(0.25,4,7);
    drawCube(0.8,0.8,0.8);
    glPopMatrix();

    //glViewport(0, 0, windowHeight, windowWidth);
    glPushMatrix();
    glTranslatef(18*0.5,0.0,18*0.5);
    glRotatef(-90,0,1,0);
    glScalef(0.25,4,7);             ///wall 01
    drawCube(0.8,0.8,0.8);
    glPopMatrix();

    glPopMatrix();

}
///...................................................................................

void display(void)
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 5000.0;
    const double a = t*90.0;

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    glFrustum(-5,5,-5,5, 4,5000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    // eye, look at, head up vector
    gluLookAt(eye_x, eye_y, eye_z, look_x, look_y, 0, 0, 1, 0);
    glViewport(0, 0, windowHeight, windowWidth);

    glPushMatrix();
for(int j=-500;j<500;j+=250){

for(int i = 0; i<50; i+=16){
    ///building_wall
        glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,ID_[0]);
    glTranslatef(0+110, i, 60+j);
    glScalef(2,2,2);
    building_wall();
    glPopMatrix();
     glDisable(GL_TEXTURE_2D);

    ///room-furniture
    glPushMatrix();
    glTranslatef(0+110, i, 60+j);
    glScalef(2,2,2);
    building_help();
    glPopMatrix();
}
}

    ///moon....
    glPushMatrix();
    glTranslatef(500,130,-1000+move_moon);
    glScalef(20,20,22);
    solid(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-1000,130,1000+move_sun);
    glScalef(25,25,25);
    solid(1,1,1);
    glPopMatrix();


    //glTranslatef(0,0,50);



    glPushMatrix();
    glTranslatef(8,5,18);
    glScalef(2,2.8,2);
    nagordola_func();
    glPopMatrix();


    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,sea_texture);
    sea();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    for(int i=-1000; i<=1000; i+=200)
    {
        glPushMatrix();


        glPushMatrix();

        glTranslatef(10.5,10,50+50+i);
        glScalef(5,5,5);
        glPushMatrix();
        glTranslatef(0,5,0);
        glRotatef(90,1,0,0);
        nagordola_stand();
        glPopMatrix();

        //other side
        glPushMatrix();
        glTranslatef(0,5,4);
        glRotatef(90,1,0,0);
        nagordola_stand2();
        glPopMatrix();

        glPopMatrix();






        glPushMatrix();


        glPushMatrix();
        glTranslatef(0,0,50-18+i);
        chorki_stand(60);
        glPopMatrix();

        glPushMatrix();
        glTranslatef(0,0,50+18+i);
        chorki_stand(60);
        glPopMatrix();

        glPopMatrix();


        glPopMatrix();
    }

    glPushMatrix();
    road();
    glPopMatrix();



    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID_[6]);
    railway_tunnel();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);



    /**
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D, ID_[3]);
        building();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);


        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D, ID_[4]);
        building2();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);



        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D, ID_[5]);
        building3();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    **/


    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D,ID_[2]);
    sky();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);



    glPushMatrix();

    glPushMatrix();
    //glColor3f(1,1,1);
    glTranslatef(31+27-3+0.5-2+0.5-2.5+6-1,13,10+4+20-4-29+5+1000+zx1);
    glRotatef(90,0,1,0);
    glScalef(0.5,0.5,0.5);
    bottleBezier(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(31+27-3+0.5-2+0.5-3,-1.5,10+4+20-4-29+1000+zx1);
    glScalef(5.5,1.8,3);
    drawCube(1,0,0);
    glPopMatrix();

    glPushMatrix();
    for(int i=0; i<=140; i+=35)
    {
        glPushMatrix();
        rail(0,i+1000+zx1,rotate1);
        glPopMatrix();
    }
    glPopMatrix();
    glPopMatrix();





    glPushMatrix();

    glPushMatrix();
    //glColor3f(1,1,1);
    glTranslatef(31+27-3+0.5-2+0.5-2.5+6-1+18,13,10+4+20-4-29+5+173.5-1000+zx2);
    glRotatef(90,0,1,0);
    glScalef(0.5,0.5,0.5);
    bottleBezier(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(31+27-3+0.5-2+0.5-3+18,-1.5,10+4+20-4-29+180-2-2-1000+zx2);
    glScalef(5.5,1.8,3);
    drawCube(1,0,0);
    glPopMatrix();

    glPushMatrix();
    for(int i=0; i<=140; i+=35)
    {
        glPushMatrix();
        rail(18,i-1000+zx2,rotate2);
        glPopMatrix();
    }
    glPopMatrix();

    glPopMatrix();









    glPushMatrix();
    railline();
    glPopMatrix();
    //rotating_bike
    glPushMatrix();
    //glRotatef( theta, axis_x, axis_y, 0.0 );
    //glRotatef( alpha, axis_x, axis_y, 0.0 );
    glTranslatef(8,5,3);
    glScalef(0.75, 0.75, 0.75);
    bike_func();
    glPopMatrix();

    //wall
    glPushMatrix();
    //glRotatef( theta, axis_x, axis_y, 0.0 );
    //glRotatef( alpha, axis_x, axis_y, 0.0 );
    deyal();
    glPopMatrix();

    glPushMatrix();
    //sky();
    glPopMatrix();

    //dheki
    for(int i=0; i<=30; i+=5)
    {
        glPushMatrix();
        //glRotatef( theta, axis_x, axis_y, 0.0 );
        //glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(15-18+i,0.2,-3+37);
        call_dheki_func();
        glPopMatrix();
    }

    //dheki
    for(int i=0; i<=40; i+=8)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(-1-5-2+i,0.2,-12);
        glRotatef(90,0,1,0);
        call_dheki_func();
        glPopMatrix();
    }

    //dolna
    for(int i=0; i<=40; i+=8)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(-8+42,0.2,13+20-i);
        glRotatef(90,0,1,0);
        call_dolna_func();
        glPopMatrix();
    }


    for(int i=0; i<=40; i+=8)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(-8-1,0.2,13+20-i);
        glRotatef(90,0,1,0);
        call_dolna_func();
        glPopMatrix();
    }

    //tree
    glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
    ///array_of_tree();
    glPopMatrix();

    //moi
    for(int i=0; i<=20; i+=5)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(2+20,0.2,-5+i);
        glRotatef(90,0,1,0);
        call_ladder();
        glPopMatrix();
    }


    //sliding moi
    for(int i=0; i<=24; i+=8)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(25,0.2,-5+i);
        ///glRotatef(90,0,1,0);
        sliding_ladder();
        glPopMatrix();
    }

    for(int i=0; i<=16; i+=8)
    {
        glPushMatrix();
//    glRotatef( theta, axis_x, axis_y, 0.0 );
//    glRotatef( alpha, axis_x, axis_y, 0.0 );
        glTranslatef(-5+2,0.2,-5+i);
        ///glRotatef(90,0,1,0);
        sliding_ladder();
        glPopMatrix();
    }

    glPopMatrix();


    /*
    //visualizing lighting effect
    glPushMatrix();
    glRotatef( theta, axis_x, axis_y, 0.0 );
    glRotatef( alpha, axis_x, axis_y, 0.0 );
    glTranslatef(-6.4, 4, 0);
    glScalef(0.25, 0.25, 0.25);
    solid(1,0,0);
    glPopMatrix();
    */


    glFlush();
    glutSwapBuffers();
}
bool stop = 1;
void all_refresh(void)
{
    AL=false;
    DL=false;
    SL=false;
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    glDisable(GL_LIGHT2);
    glDisable(GL_LIGHT3);
    glDisable(GL_LIGHT4);
    glDisable(GL_LIGHT5);
}
void refresh(void)
{

    glDisable(GL_LIGHT2);
    glDisable(GL_LIGHT3);
    glDisable(GL_LIGHT4);
    glDisable(GL_LIGHT5);
}
void sokal(void)
{
    refresh();
    AL= true;
    DL= true;
    SL= true;
    glEnable(GL_LIGHT2);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    move_sun = 0;
    move_moon = -3000;
    light();

}
void dupur(void)
{
    refresh();
    AL= true;
    DL= true;
    SL= true;
    glEnable(GL_LIGHT3);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    move_sun = -1000;
    move_moon = -3000;
    light();

}
void bikal(void)
{
    refresh();
    AL= true;
    DL= true;
    SL= true;
    glEnable(GL_LIGHT4);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    move_sun = -2000;
    move_moon = -3000;
    light();

}
void raat(void)
{
    refresh();
    AL= true;
    DL= true;
    SL= true;
    glEnable(GL_LIGHT5);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    move_sun = -3000;
    move_moon = 0;
    light();

}
void myKeyboardFunc( unsigned char key, int x, int y )
{
    switch ( key )
    {


    case '*':
        stop = 1-stop;

        break;
    case '0':
        fanRotate = !fanRotate;
        uRotate = false;
        //bRotate = false;
        f_x=0.0;
        f_y=1.0;
        break;
    //movement (front, back, left, right, up, down)
    case 'w':
    case 'W':
        eye_z-=0.5;
//        if(eye_z<1)
//            eye_z = 10;
        break;
    case 's':
    case 'S':
        eye_z+=0.5;
//        if(eye_z>20)
//            eye_z = 10;
        break;
    case 'a':
    case 'A':
        eye_x-=0.5;
//        if(eye_x<-12)
//            eye_x = 2;
        break;

    case 'd':
    case 'D':
        eye_x+=0.5;
//        if(eye_x>10)
//            eye_x = 2;
        break;

    case 'u':
    case 'U':
        eye_y+=0.5;
//        if(eye_y>8)
//            eye_y = 2;
        break;

    case 'p':
    case 'P':
        eye_y-=0.5;
//        if(eye_y<-2)
//            eye_y = 2;
        break;

    //look(up, down, left, right)
    case 'i':
    case 'I':
        look_y+=0.8;
//        if(look_y>5)
//            look_y = 0;
        break;
    case 'k':
    case 'K':
        look_y-=0.8;
//        if(look_y<-3)
//            look_y = 0;
        break;

    case 'j':
    case 'J':
        look_x-=0.8;
//        if(look_x<-20)
//            look_x = 2;
        break;

    case 'l':
    case 'L':
        look_x+=0.8;
//        if(look_x>20)
//            look_x = 2;
        break;

    case 'h':
    case 'H':
        bRotate = !bRotate;
        uRotate = false;
        axis_x=0.0;
        axis_y=1.0;
        break;

    case 'f':
    case 'F':
        dRotate = !dRotate;
        uRotate = false;
        bRotate = false;
        d_x = 1.0;
        d_y = 0.0;
        break;


    case 'v':
    case 'V':
        uRotate = !uRotate;
        bRotate = false;
        axis_x=1.0;
        axis_y=0.0;
        break;

    case '1':
        AL = !AL;
        light();
        break;
    case '2':
        DL = !DL;
        light();
        break;
    case '3':
        SL = !SL;
        light();
        break;

    case '4':
        L0 = !L0;
        if(L0)
            glEnable(GL_LIGHT0);
        else
            glDisable(GL_LIGHT0);
        break;

    case '5':
        L1 = !L1;
        if(L1)
            glEnable(GL_LIGHT1);
        else
            glDisable(GL_LIGHT1);
        break;

    case '6':
        L2 = !L2;
        if(L2)
        {
            glEnable(GL_LIGHT2);

        }
        else
        {

            glDisable(GL_LIGHT2);
        }
        break;
    case '7':
        L3 = !L3;
        if(L3)
        {
            glEnable(GL_LIGHT3);

        }
        else
        {
            glDisable(GL_LIGHT3);
        }
        break;
    case '8':
        L4 = !L4;
        if(L4)
        {
            glEnable(GL_LIGHT4);

        }
        else
        {
            glDisable(GL_LIGHT4);
        }
        break;
    case '9':
        L5 = !L5;
        if(L5)
        {
            glEnable(GL_LIGHT5);

        }
        else
        {
            glDisable(GL_LIGHT5);
        }
        break;
    case '?':
        all_refresh();
        break;
    case '!':
        sokal();
        break;
    case '@':
        dupur();
        break;
    case '#':
        bikal();
        break;
    case '$':
        raat();
        break;

    case 27:	// Escape key
        exit(1);
    }
}
void animate()
{

    light();


    if(sea_texture==ID_[9])
        sea_texture = ID_[10];
    else if(sea_texture == ID_[10])
        sea_texture = ID_[11];
    else if(sea_texture == ID_[11])
        sea_texture = ID_[12];
    else
        sea_texture = ID_[9];


    rotate1+=2;
    rotate2-=2;
    if(rotate1>360)
        rotate1=0;
    if(rotate2<-360)
        rotate2=0;

    if(stop)
    {
        zx1-=1;
        zx2+=1;
        if(zx1<-2000)
            zx1=0.0;
        if(zx2>2000)
            zx2=0.0;
    }
    if(dRotate== true)
    {
        static bool ch = false, ch2 = true;
        if(rot<45.0 && ch2)
        {
            rot += 2;
        }
        else
            ch = true;
        if(ch && rot>-45.0)
        {
            rot -= 2;
            ch2 = false;
        }
        else
        {
            ch2 = true;
            ch = false;
        }
    }




    if (bRotate == true)
    {
        theta += 0.1;
        if(theta > 360.0)
            theta -= 360.0*floor(theta/360.0);
    }

    if (uRotate == true)
    {
        alpha += 0.1;
        if(alpha > 360.0)
            alpha -= 360.0*floor(alpha/360.0);
    }
    if (fanRotate == true)
    {
        fan_rot -= 2;
        if(fan_rot <= 0.0)
            fan_rot  = 360.0;
    }

    angle+=0.50;
    if(angle>360)
        angle-=360;

    glutPostRedisplay();

}

void LoadTexture(const char*filename, unsigned int* ID)
{


    glGenTextures(1, ID);
    glBindTexture(GL_TEXTURE_2D, *ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, *ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

int main (int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowPosition(100,100);
    glutInitWindowSize(windowHeight, windowWidth);
    glutCreateWindow("Kids Park Project");
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\brick.bmp", &ID_[0]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\grass.bmp", &ID_[1]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\sky.bmp", &ID_[2]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\building.bmp", &ID_[3]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\building2.bmp", &ID_[4]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\building3.bmp", &ID_[5]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\railway_tunnel.bmp", &ID_[6]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\road.bmp", &ID_[7]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\road1.bmp", &ID_[8]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\sea1.bmp", &ID_[9]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\sea2.bmp", &ID_[10]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\sea3.bmp", &ID_[11]);
    LoadTexture("C:\\Users\\kazik\\OneDrive\\Desktop\\project4_2\\sea4.bmp", &ID_[12]);


    sea_texture = ID_[9];
    /***
    printf("Enter below key for controlling....");
    printf("\n h/H horizontal rotation\n v/V vertical rotation");
    printf("\n 1  light1(SpotLight) RED\n 2 light2(LEFT) GREEN\n 3 light3(RIGHT) BLUE");
    printf("\n 4  ambient\n 5 diffuse\n 6 specular\n");

    printf("\nEnter below key for controlling movement....fixed");
    printf("\n w/W  front\t s/S back");
    printf("\n a/A  left\t d/D right");
    printf("\n u/U  up\t p/P down");

    printf("\n\nEnter below key for controlling look....");
    printf("\n i/I up\t\t k/K down");
    printf("\n j/J  left\t l/L right");
    printf("\n\nEnter Esc for quit\n");
       **/

    glShadeModel( GL_SMOOTH );
    glEnable( GL_DEPTH_TEST );
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);



    cout<<"AL: "<<AL<<" "<<"DL: "<<DL<<"SL: :"<<SL<<endl;
    cout<<"Light2: "<<L2<<" "<<"Light3: "<<L3<<endl;

    glutKeyboardFunc(myKeyboardFunc);
    glutDisplayFunc(display);
    glutIdleFunc(animate);
    glutMainLoop();

    return 0;
}





