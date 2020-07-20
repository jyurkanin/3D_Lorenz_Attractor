#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <deque>
#ifdef Success //this is stupid
  #undef Success
#endif
#include <eigen3/Eigen/Dense>


#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080
#define STEREO_DISPARITY .068 //milimeters.
#define FOCAL_LEN .05 //5cm idk bruh

#define EPS 1e-7

using namespace Eigen;

Display *dpy;
Window w;
GC gc;

pthread_t w_thread;
int is_window_alive;

typedef struct {
    int xl, yl, xr, yr;
    float theta_l, theta_r, phi_l, phi_r; //not sure if angular or rectangular coords are better for displaying 
} ProjectedPoints;

typedef struct{
    int xl, yl, xr, yr;
} StereoPixels;

typedef struct{
    float x, y, z;
} Point3D;

typedef struct{
    Vector3f pos;
    float radius;
} Sphere;

typedef struct{
    Vector3f pos;
    Matrix3f rot;
    Vector3f vertices;
} RigidBody;

Vector3f origin = {0,0,0};

Matrix3f get_rotation(float x, float y, float z){
    Matrix3f rotx; rotx << 1,0,0, 0,cosf(x),-sinf(x), 0,sinf(x),cosf(x);
    Matrix3f roty; roty << cosf(y),0,sinf(y), 0,1,0, -sinf(y),0,cosf(y);
    Matrix3f rotz; rotz << cosf(z),-sinf(z),0, sinf(z),cosf(z),0, 0,0,1;

    return rotx*roty*rotz;
}

ProjectedPoints project_point(Vector3f point_v, Vector3f view_v){
    Point3D view;
    Point3D point;
    ProjectedPoints pp;

    view.x = view_v[0];
    view.y = view_v[1];
    view.z = view_v[2];
    
    point.x = point_v[0];
    point.y = point_v[1];
    point.z = point_v[2];
    
    pp.theta_l = atan2f(point.y - view.y, point.x - view.x);
    pp.theta_r = atan2f(point.y - view.y - STEREO_DISPARITY, point.x - view.x);
    pp.phi_l = atan2f(point.z - view.z, point.x - view.x);
    pp.phi_r = atan2f(point.z - view.z, point.x - view.x);

    /*
    pp.xl = FOCAL_LEN*(point.x - view.x)/(point.y - view.y);
    pp.yl = FOCAL_LEN*(point.z - view.z)/(point.y - view.y);
    pp.xr = FOCAL_LEN*(point.x - view.x - STEREO_DISPARITY)/(point.y - view.y);
    pp.yr = FOCAL_LEN*(point.z - view.z)/(point.y - view.y);
    */
    return pp;
}

StereoPixels projection_to_pixels(ProjectedPoints pp){
    StereoPixels sp;
    float theta_max = M_PI*.7;
    float phi_max = M_PI*.7;
    float scale_x = (SCREEN_WIDTH/4)/theta_max;
    float scale_y = (SCREEN_HEIGHT/2)/phi_max;

//    printf("angle %f %f\n", pp.theta_l, pp.phi_l);
    if(abs(pp.theta_l) > (theta_max + EPS) || abs(pp.phi_l) > (phi_max + EPS)){
        sp.xl = -1;
        sp.yl = -1;
    }
    else{
        sp.xl = floor(pp.theta_l*scale_x + SCREEN_WIDTH/4);
        sp.yl = floor(-pp.phi_l*scale_y + SCREEN_HEIGHT/2);
    }
    
    if(abs(pp.theta_r) > (theta_max + EPS) || abs(pp.phi_r) > (phi_max + EPS)){
        sp.xr = -1;
        sp.yr = -1;
    }
    else{
        sp.xr = floor(pp.theta_r*scale_x + 3*SCREEN_WIDTH/4);
        sp.yr = floor(-pp.phi_r*scale_y + SCREEN_HEIGHT/2);
    }
    return sp;
}

void draw_stereo_line(Vector3f start_point, Vector3f end_point){
    //should be a line sloping to the top right.

    ProjectedPoints pp_start = project_point(start_point, origin);
    ProjectedPoints pp_end = project_point(end_point, origin);
    
    StereoPixels sp_start = projection_to_pixels(pp_start);
    StereoPixels sp_end = projection_to_pixels(pp_end);

    if(!(sp_start.xl == -1 || sp_end.xl == -1))
        XDrawLine(dpy, w, gc, sp_start.xl, sp_start.yl, sp_end.xl, sp_end.yl);
    if(!(sp_start.xr == -1 || sp_end.xr == -1))
        XDrawLine(dpy, w, gc, sp_start.xr, sp_start.yr, sp_end.xr, sp_end.yr);
}

void draw_stereo_point(Vector3f point){
    ProjectedPoints pp = project_point(point, origin);
    StereoPixels sp = projection_to_pixels(pp);
    XDrawPoint(dpy, w, gc, sp.xl, sp.yl);
    XDrawPoint(dpy, w, gc, sp.xr, sp.yr);
}

void draw_stereo_sphere(Sphere sphere){
    ProjectedPoints pp = project_point(sphere.pos, origin);
    StereoPixels sp = projection_to_pixels(pp);
    Vector3f temp(sphere.pos);
    float radius_angle = atan2f(sphere.radius, temp.norm())*128;
    XFillArc(dpy, w, gc, sp.xl, sp.yl, radius_angle, radius_angle, 0, 360*64);
    XFillArc(dpy, w, gc, sp.xr, sp.yr, radius_angle, radius_angle, 0, 360*64);
}

int rainbow(int c){
  static int red = 0x00;
  static int green = 0x00;
  static int blue = 0xFF;
  static int state = 0;
  static int counter = 0;
  
  counter++;
  if(counter < c){
      int out = 0;
      out = (red << 16) | (blue << 8) | green;
      return out;
  }
  
  counter = 0;

  switch(state){
  case 0:
      green++;
      if(green == 0x8F) state = 1;
      break;
  case 1:
      blue--;
      if(blue == 0x00) state = 2;
      break;
  case 2:
      blue++;
      if(blue == 0x8F) state = 3;
      break;
  case 3:
      green--;
      if(green == 0x00) state = 0;
      break;
  }

  int out = 0;
  out = (red << 16) | (green << 8) | blue;
  return out;
  
}

void lorenz_ode(float *X, float *Xd){
    float sigma = 10;
    float beta = 8/3;
    float row = 22;
    
    Xd[0] = sigma*(X[1]-X[0]);
    Xd[1] = X[0]*(row-X[2])-X[1];
    Xd[2] = X[0]*X[1]-beta*X[2];
}

void runge_kutta(float *X, float *Xt1, void (ode)(float*,float*), float ts, int len){
    float temp[3];
    float k1[3]; ode(X, k1); //Calculate derivative.
    
    for(int i = 0; i < len; i++) temp[i] = X[i]+.5*ts*k1[i];
    float k2[3]; ode(temp, k2);
    
    for(int i = 0; i < len; i++) temp[i] = X[i]+.5*ts*k2[i];
    float k3[3]; ode(temp, k3);
    
    for(int i = 0; i < len; i++) temp[i] = X[i]+ts*k3[i];
    float k4[3]; ode(temp, k4);
    
    for(int i = 0; i < len; i++){
        Xt1[i] = X[i] + (ts/6)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

void* window_thread(void*){
    XSetForeground(dpy, gc, 0xFF);

    float num = .07;
    
    Vector3f d(3,0,0);
    
    float X[3] = {10,.1,.2};
    float new_X[3];

    Matrix3f rot = get_rotation(0,0,0);
    Matrix3f i_rot = get_rotation(0,0,.002);
    
    std::deque<Vector3f> lines;     

    for(int j = 0; j < 10000; j++){
//        d[0] += .08;
        rot = rot*i_rot;
        XSetForeground(dpy, gc, 0);
        XFillRectangle(dpy, w, gc, 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
        for(int i = 0; i < 10; i++){
            runge_kutta(X, new_X, lorenz_ode, .01, 3);
            lines.push_back(Vector3f(new_X)); //create new vector3f from array
            memcpy(X, new_X, 3*sizeof(float));
            
            if(lines.size() > 300){
                lines.pop_front();
            }
            
            XSetForeground(dpy, gc, rainbow(8));
            for(int k = 1; k < lines.size(); k++){
                draw_stereo_line((rot*lines[k-1]*num)+d, (rot*lines[k]*num)+d);
            }
            XFlush(dpy);
//            usleep(10);
        }
    }
    return 0;
}


void init_window(){
    dpy = XOpenDisplay(":0.0");
    w = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, 0, 0);

    Atom wm_state   = XInternAtom (dpy, "_NET_WM_STATE", true );
    Atom wm_fullscreen = XInternAtom (dpy, "_NET_WM_STATE_FULLSCREEN", true );
    XChangeProperty(dpy, w, wm_state, XA_ATOM, 32, PropModeReplace, (unsigned char *)&wm_fullscreen, 1);
    
    XSelectInput(dpy, w, StructureNotifyMask | ExposureMask | KeyPressMask);
    XClearWindow(dpy, w);
    XMapWindow(dpy, w);
    gc = XCreateGC(dpy, w, 0, 0);
    is_window_alive = 1;
    
    XEvent e;
    do{
        XNextEvent(dpy, &e);        
    } while(e.type != MapNotify);

    pthread_create(&w_thread, NULL, &window_thread, NULL);    
}



void del_window(){
    XDestroyWindow(dpy, w);
    XCloseDisplay(dpy);
    is_window_alive = 0;
    return;
}



int main(){
    init_window();
    while(is_window_alive){
        
    }
    del_window();
    
}
