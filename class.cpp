//0: xyz grid reference
//in Km, and nanoseconds
using namespace std;

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <fstream>

long double C = 300000000.0;
long double R = 6400.0;
long double TIMERESOLUTION = 100.0;
// in nanoseconds, simulation-time for a tick
long double RTIMERESOLUTION = 100000.0;
// in microseconds, real-time for a tick
#define ERRORRATE 1000
#define ALTCOEFF 200
#define ERRORTOLORATE 10000
#define PRECISIONDEFINITION cout.precision(15)
#define PRINTLINE cout << "---------------------------------" << endl

//string FILENAME = "log";
//ofstream fout( FILENAME.c_str(), std::ios::app);

//class sate;
//class center;
//class user;

class sate;
class center;

vector< sate > vSate;
vector< center > vCenter;

long double tpr_x(long double theta, long double phi, long double rad){
  return rad*sin(theta)*cos(phi);
}
long double tpr_y(long double theta, long double phi, long double rad){
  return rad*sin(theta)*sin(phi);
}
long double tpr_z(long double theta, long double phi, long double rad){
  return rad*cos(theta);
}
long double randomFloat(long double a, long double b) {
    long double random = ((long double) rand()) / (long double) RAND_MAX;
    long double diff = b - a;
    long double r = random * diff;
    return a + r;
}
long double dist(long double x, long double y, long double z, long double a, long double b, long double c){
  return sqrt((x-a)*(x-a)+(y-b)*(y-b)+(z-c)*(z-c));
}
long double cDist(long double x, long double y, long double z, long double a, long double b, long double c, long long error){
  return dist(x,y,z,a,b,c)-(((long double)(C * error)+0.0)/1000000000);
}
int chkComm(long double x,long double y,long double z,long double a,long double b, long double c){
  long double check=(a*x-x*x+b*y-y*y+c*z-z*z)
                    /(sqrt((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))*sqrt(x*x+y*y+z*z));
//  cout << "DEBUG: COS(THETA) = "<< check << endl;
//  cout << "DEBUG: A = "<< (a*x-x*x+b*y-y*y+c*z-z*z) << endl;
//  cout << "DEBUG: /B = "<< (sqrt((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))*sqrt(x*x+y*y+z*z)) << endl;
//  cout << "DEBUG: /p = "<< sqrt((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z)) << endl;
//  cout << "DEBUG: /q = "<< sqrt(x*x+y*y+z*z) << endl;
//  cout << "DEBUG: A/B = "<< (a*x-x*x+b*y-y*y+c*z-z*z)/(sqrt((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))*sqrt(x*x+y*y+z*z)) << endl;
  if(0.0<=check){
    return 1;
  }
  else{
    return 0;
  }
}
class xyz{
public:
  long double x_;
  long double y_;
  long double z_;

  xyz(){}
};

class sate{
private:
  long double x_;
  long double y_;
  long double z_;

  long double theta_;
  long double phi_;
  long double rad_;

  long long error_;

public:
  sate(){}
  void setSateposition(){
    theta_ = randomFloat(0.0, 2*M_PI);
    phi_ =  randomFloat(0.0, M_PI);
    rad_ = R+randomFloat(100,100+ALTCOEFF);

    x_=tpr_x(theta_, phi_, rad_);
    y_=tpr_y(theta_, phi_, rad_);
    z_=tpr_z(theta_, phi_, rad_);
    cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
  }
/*  ON DEV
    void setSateposition(int i, long double k){
    if(i > (vCenter.size()-1)){
      cout << "FATAL ERROR: THERE ARE UP TO #" << vCenter.size()-1 << " CENTERS" << endl;
    }
    else{
      center* Center = vCenter[i];

    }
  }
  */
  void setSateposition(int i, long double p, long double q);
  void setSateposition(long double a, long double b, long double c){
      theta_ = a;
      phi_ = b;
      rad_ = c;

      x_=tpr_x(theta_, phi_, rad_);
      y_=tpr_y(theta_, phi_, rad_);
      z_=tpr_z(theta_, phi_, rad_);
      cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
    }
  void setError(){
    error_=error_+TIMERESOLUTION*ERRORRATE*(rand()%100);
  }
  void setError(long long error){
    error_=error;
  }
  long double getx(void){
    return x_;
  }
  long double gety(void){
    return y_;
  }
  long double getz(void){
    return z_;
  }
  long double gett(void){
    return theta_;
  }
  long double getp(void){
    return phi_;
  }
  long double getr(void){
    return rad_;
  }
  long long getError(void){
    return error_;
  }

};
class center{
private:
  long double x_;
  long double y_;
  long double z_;

  long double theta_;
  long double phi_;
  long double rad_;

  long long tError_;

  vector <sate*> sate_;

public:
  center(long long tError){
    tError_ = tError;
  }
  void setCenterposition(long double a, long double b, long double c){
      theta_ = a;
      phi_ = b;
      rad_ = c;

      x_=tpr_x(theta_, phi_, rad_);
      y_=tpr_y(theta_, phi_, rad_);
      z_=tpr_z(theta_, phi_, rad_);
      cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
    }
  void setCenterposition(){
    theta_ = randomFloat(0.0, 2*M_PI);
    phi_ =  randomFloat(0.0, M_PI);
    rad_ = R;

    x_=tpr_x(theta_, phi_, rad_);
    y_=tpr_y(theta_, phi_, rad_);
    z_=tpr_z(theta_, phi_, rad_);
    cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
  }
  void init(){
    cout << "DEBUG: INITIALIZING CENTER WITH (NEW) SATES SET" << endl;
    sate_.clear();
    for(int i=0; i < vSate.size(); i++){
      //if coummunicatable,
      sate* Sate = &vSate[i];
      if(chkComm(x_,y_,z_,Sate->getx(),Sate->gety(),Sate->getz())){
        sate_.push_back(Sate);
        cout << "DEBUG: A COMMUNCATABLE SATE FOUND WITH XYZ: "<<Sate->getx()<<" "<<Sate->gety()<<" "<<Sate->getz()<<endl;
      }
      else{
        cout << "DEBUG: AN UNCOMMUNICATABLE SATE FOUND WITH XYZ: "<<Sate->getx()<<" "<<Sate->gety()<<" "<<Sate->getz()<<endl;
      }
    }
  }
  void syncTime(){
    for(int i=0; i < sate_.size(); i++){
      sate* Sate = sate_[i];
      if(abs(Sate->getError())>tError_){
        Sate->setError(.0);
      }
    }
  }
  long double getx(void){
    return x_;
  }
  long double gety(void){
    return y_;
  }
  long double getz(void){
    return z_;
  }
  long double gett(void){
    return theta_;
  }
  long double getp(void){
    return phi_;
  }
  long double getr(void){
    return rad_;
  }
};

class user{
private:
  long double x_;
  long double y_;
  long double z_;

  long double theta_;
  long double phi_;
  long double rad_;

  vector<sate*> commSate_;

public:
  user(){}
  void init(){
    for(int i=0;i<vSate.size();i++){
      sate* Sate = &vSate[i];
      if(chkComm(x_,y_,z_,Sate->getx(),Sate->gety(),Sate->getz())){
        commSate_.push_back(Sate);
        cout << "DEBUG: COMMUNICATABLE SATE FOUUND WITH TPR: "<<Sate->getx()<<" "<<Sate->gety()<<" "<<Sate->getz()<<endl;
      }
    }
  }
  void setUserposition(long double a, long double b, long double c){
    theta_ = a;
    phi_ = b;
    rad_ = c;

    x_=tpr_x(theta_, phi_, rad_);
    y_=tpr_y(theta_, phi_, rad_);
    z_=tpr_z(theta_, phi_, rad_);
    cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
  }
  void setUserposition(){
    theta_ = randomFloat(0.0, 2*M_PI);
    phi_ =  randomFloat(0.0, M_PI);
    rad_ = R;

    x_=tpr_x(theta_, phi_, rad_);
    y_=tpr_y(theta_, phi_, rad_);
    z_=tpr_z(theta_, phi_, rad_);
    cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<" XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
  }

  void printCoordinate(void){
    if(commSate_.size() < 3){
      cout << "FATAL ERROR: NOT ENOUGH SATES" << endl;
      cout << "DEBUG: YOU HAVE "<<commSate_.size()<<" COMMUNICATABLE SATES NOW" << endl;
    }
    else{
      int i, j, k, n;
      long double ax, ay, az, at, ap, ar;
      ax = ay = az = at = ap = ar = 0.0;
      for(i=0, n=0;i<commSate_.size()-2;i++){
        for(j=i+1;j<commSate_.size()-1;j++){
          for(k=j+1;k<commSate_.size();k++){
            // getCoordinate();
            // vector[i]와 j, k를 이용하여 getCoordinate한 것들을, 평균
              sate* S1 = commSate_[i];
              sate* S2 = commSate_[j];
              sate* S3 = commSate_[k];
              cout << "NOW CALCULATING WITH "<<i<<", "<<j<<", "<<k<<" SATES" << endl;
              xyz* XYZ = getCoordinate(S1->getx(),S1->gety(),S1->getz(),
                                       S2->getx(),S2->gety(),S2->getz(),
                                       S3->getx(),S3->gety(),S3->getz(),
                                       cDist(x_,y_,z_,S1->getx(),S1->gety(),S1->getz(),S1->getError()),
                                       cDist(x_,y_,z_,S2->getx(),S2->gety(),S2->getz(),S2->getError()),
                                       cDist(x_,y_,z_,S3->getx(),S3->gety(),S3->getz(),S3->getError()));
            //A의 좌표, B의 좌표, C의 좌표, A까지 추정 거리, B까지 추정 거리, C까지 추정 거리
            cout << XYZ->x_<<" "<< XYZ->y_<<" "<< XYZ->z_<<" "<< endl;
            ax = ax + XYZ->x_;
            ay = ay + XYZ->y_;
            az = az + XYZ->z_;

            delete XYZ;
            n++;
          }
        }
      }

      ax = ax / (n+.0);
      ay = ay / (n+.0);
      az = az / (n+.0);

      PRINTLINE;
      cout << "AVERAGE X ASSUMPTION: " << ax << endl;
      cout << "AVERAGE Y ASSUMPTION: " << ay << endl;
      cout << "AVERAGE Z ASSUMPTION: " << az << endl;
      PRINTLINE;
      at = acos(az/sqrt(ax*ax+ay*ay+az*az));
      ap = atan2(ay,ax);
      ar = sqrt(ax*ax+ay*ay+az*az);
      cout << "AVERAGE THETA ASSUMPTION: " << at << endl;
      cout << "AVERAGE PHI ASSUMPTION: " << ap << endl;
      cout << "AVERAGE RHO ASSUMPTION: " << ar << endl;
    }
  }

  xyz* getCoordinate(long double a, long double b, long double c, long double p, long double q, long double r, long double i, long double j, long double k, long double d, long double s, long double l){

    long double x_, y_, z_;
    x_ = -(((4*c*q - 4*b*r)*(2*k*(a*a + b*b + c*c - d*d + R*R) -
               2*c*(i*i + j*j + k*k - l*l + R*R)) - (4*c*j -
               4*b*k)*(2*r*(a*a + b*b + c*c - d*d + R*R) -
               2*c*(p*p + q*q + r*r + R*R - s*s)))/((-16*(c*c))*j*p +
            16*b*c*k*p + 16*(c*c)*i*q - 16*a*c*k*q - 16*b*c*i*r +
            16*a*c*j*r));

    y_ = -((c*i*i*p + c*j*j*p - a*a*k*p - b*b*k*p - c*c*k*p + d*d*k*p +
          c*k*k*p - c*l*l*p - c*i*p*p + a*k*p*p - c*i*q*q + a*k*q*q +
          a*a*i*r + b*b*i*r + c*c*i*r - d*d*i*r - a*i*i*r - a*j*j*r -
          a*k*k*r + a*l*l*r - c*i*r*r + a*k*r*r - c*i*R*R + a*k*R*R +
          c*p*R*R - k*p*R*R - a*r*R*R + i*r*R*R + c*i*s*s -
          a*k*s*s)/(2*(-c*j*p + b*k*p + c*i*q - a*k*q - b*i*r +
            a*j*r)));

    z_ = (-b*i*i*p + a*a*j*p + b*b*j*p + c*c*j*p -d*d*j*p - b*j*j*p -
        b *k*k* p + b*l*l*p + b*i*p*p - a*j*p*p - a*a*i*q - b*b*i*q -
        c*c* i* q + d*d*i*q + a*i*i*q + a*j*j*q + a*k*k*q - a*l*l*q +
        b*i*q*q - a*j*q*q + b*i*r*r - a*j*r*r + b*i*R*R - a*j*R*R -
        b*p*R*R + j*p*R*R + a*q*R*R - i*q*R*R - b*i*s*s +
        a*j*s*s)/(2*(c*j*p - b*k*p - c*i*q + a*k*q + b*i*r - a*j*r));
    xyz* XYZ = new xyz();
    XYZ->x_=x_;
    XYZ->y_=y_;
    XYZ->z_=z_;

    return XYZ;
  }
  long double getx(void){
    return x_;
  }
  long double gety(void){
    return y_;
  }
  long double getz(void){
    return z_;
  }
  long double gett(void){
    return theta_;
  }
  long double getp(void){
    return phi_;
  }
  long double getr(void){
    return rad_;
  }
};

/* REFERENCE TO CHECK WHETHER IMPLEMENTATION WAS CORRECTLY DONE

Solve[{(x - a)^2 + (y - b)^2 + (z - c)^2 - (x^2 + y^2 + z^2) ==
   d^2 - R^2, (x - p)^2 + (y - q)^2 + (z - r)^2 - (x^2 + y^2 + z^2) ==
    s^2 - R^2, (x - i)^2 + (y - j)^2 + (z - k)^2 - (x^2 + y^2 +
      z^2) == l^2 - R^2}, {x, y, z}]

{{x -> -(((4 c q - 4 b r) (2 k (a^2 + b^2 + c^2 - d^2 + R^2) -
           2 c (i^2 + j^2 + k^2 - l^2 + R^2)) - (4 c j -
           4 b k) (2 r (a^2 + b^2 + c^2 - d^2 + R^2) -
           2 c (p^2 + q^2 + r^2 + R^2 - s^2)))/(-16 c^2 j p +
        16 b c k p + 16 c^2 i q - 16 a c k q - 16 b c i r +
        16 a c j r)),
  y -> -((c i^2 p + c j^2 p - a^2 k p - b^2 k p - c^2 k p + d^2 k p +
        c k^2 p - c l^2 p - c i p^2 + a k p^2 - c i q^2 + a k q^2 +
        a^2 i r + b^2 i r + c^2 i r - d^2 i r - a i^2 r - a j^2 r -
        a k^2 r + a l^2 r - c i r^2 + a k r^2 - c i R^2 + a k R^2 +
        c p R^2 - k p R^2 - a r R^2 + i r R^2 + c i s^2 -
        a k s^2)/(2 (-c j p + b k p + c i q - a k q - b i r +
          a j r))),
  z -> (-b i^2 p + a^2 j p + b^2 j p + c^2 j p - d^2 j p - b j^2 p -
      b k^2 p + b l^2 p + b i p^2 - a j p^2 - a^2 i q - b^2 i q -
      c^2 i q + d^2 i q + a i^2 q + a j^2 q + a k^2 q - a l^2 q +
      b i q^2 - a j q^2 + b i r^2 - a j r^2 + b i R^2 - a j R^2 -
      b p R^2 + j p R^2 + a q R^2 - i q R^2 - b i s^2 +
      a j s^2)/(2 (c j p - b k p - c i q + a k q + b i r - a j r))}}



coordinate      distance
(a,b,c)        d
(p,q,r)        s
(i,j,k)        l

*/

void sate::setSateposition(int i, long double p, long double q){
  if(i > (vCenter.size()-1)){
    cout << "FATAL ERROR: THERE ARE UP TO #" << vCenter.size()-1 << " CENTERS" << endl;
  }
  else{
    long double st, sp, sr, sx, sy, sz;
    center* Center = &vCenter[i];

    while(1){
      sr = randomFloat(p, q);
      st = randomFloat(0, 2*M_PI);
      sp = randomFloat(0, M_PI);

      sx = tpr_x(st, sp, sr);
      sy = tpr_y(st, sp, sr);
      sz = tpr_z(st, sp, sr);

//    cout << "DEBUG: A TEMP SATE GENERATED WITH TPR: "<<st<<" "<<sp<<" "<<sr<<endl;
      if(chkComm(Center->getx(),Center->gety(),Center->getz(),sx,sy,sz)){
        rad_=sr;
        theta_=st;
        phi_=sp;

        x_=sx;
        y_=sy;
        z_=sz;

        cout << "DEBUG: A SATE GENERATED  WITH ";
        break;
      }
      else{

      }
    }
    cout << "TPR: "<<theta_<<" "<<phi_<<" "<<rad_<<endl;
    cout << "XYZ: "<<x_<<" "<<y_<<" "<<z_<<endl;
  }
}


int main(){

  PRECISIONDEFINITION;
  srand(time(NULL));
  int i;
  long double a, b, c;

  system("clear");
  cout << "THIS IS A TEST PROGRAMME" << endl;
  cout << "SETTINGS: "<< endl << "TIME IN SIMULATION FOR A TICK(n): " << TIMERESOLUTION << endl;
  cout << "TIME IN REAL WORLD FOR A TICK(u): " << RTIMERESOLUTION << endl;
/*
  cout << "HOW MANY SATES ARE THERE IN THIS SYSTEM: " << endl;
  cin >> i;

  for(int k=0;k<i;k++){
    vSate.push_back(sate());
  }

//  cout << Sate.getx()<<" "<<" "<<Sate.gety()<<" "<<Sate.getz()<<endl;

  cout << vSate.size() << " OF SATES HAVE BEEN CREATED" << endl;
*/
  center* Center;
  sate* Sate;

  cout << "HOW MANY CENTERS ARE THERE IN THIS SYSTEM: " << endl;
  cin >> i;
  if(i<1){
    system("clear");
    return 0;
  }

  for(int k=0;k<i;k++){
    PRINTLINE;
    vCenter.push_back(center(ERRORTOLORATE));
    Center = &vCenter[k];
    cout << "INPUT CENTER "<<k<<"'S POSITION IN TPR <THETA, PHI>: "<<endl;
//    cin >> a >> b >> c;
    cin >> a >> b;
    Center->setCenterposition(a,b,R);
    cout <<"CENTER "<<k<<" GENERATED" << endl;
  }
  for(int k=0;;k++){
    PRINTLINE;
/*
    cout << "INPUT "<<k<<" SATE'S POSITION VALUE: "<<endl;
    vSate.push_back(sate());
    Sate = &vSate[k];
    cin >> a >> b >> c;
    if(a < 0 || b < 0 || c < 0){
      break;
    }
    Sate->setSateposition(a,b,c);
*/
    cout << "INPUT SATE "<<k<<"'S TARGET CENTER, MIN RAD, MAX RAD: "<<endl;
    cout << "TO SKIP, INPUT ANY NEGATIVE INTEGERS" << endl;
    vSate.push_back(sate());
    Sate = &vSate[k];

    int n;
    cin >> n >> a >> b;
    if(n < 0 || b < 0 || c < 0){
      break;
    }
    Sate->setSateposition(n, a, b);
    Sate->setError(0.0);

    for(int j=0;j<vCenter.size();j++){
      cout << "CENTER #"<<j<<" SELECTED"<<endl;
      Center = &vCenter[j];
      Center->init();
    }
  }
/*
  for(int k=0;k<i;k++){
    Sate = &vSate[i];
    cin >> a >> b >> c;
    Sate->setSateposition(a, b, c);
    Sate->setError(.0);
  }
*/
  user* User = new user();
  PRINTLINE;
  cout << "INPUT USER'S POS IN TPR <THETA, PHI>: ";
//  cin >> a >> b >> c;
  cin >> a >> b;
  User->setUserposition(a, b ,R);
  User->init();
  User->printCoordinate();

/*
  //Time Simualtion
  long double elapsed = .0;
  while(1){
    cout << "TIME OF SIMULATION ELAPSED IN SECONDS: "<<(long double)elapsed<<endl;
    elapsed=elapsed+(long double)(0.0+TIMERESOLUTION/1000000000);
    usleep(RTIMERESOLUTION);
  }

 */
}
