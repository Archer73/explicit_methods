#include <stdio.h>
#include "adams.h"
#include "adams5.h"
#include "rk4.h"
#include "rk5.h"

#define EXIT_IF_0(X) if(!(X)) goto error

#define EQUATIONS_NUM 2

#define STEP 1.E-3

enum VAL_NAME {V, X};

struct user_data{
  double k; //spring const
  double m; //pendulum mass
};

double RightSideV(double const x, double const *y, void *userdata){
  struct user_data *data = userdata;
  return -data->k / data->m * y[X]; //Hooke's law
}

double RightSideX(double const x, double const *y, void *userdata){
  return y[V];
}

int TestAdams(void){
  FILE * a_res = fopen("adams.txt", "w");
  a_data *adams_data;
  EXIT_IF_0(AdamsInitData(&adams_data, EQUATIONS_NUM));
  double vals[EQUATIONS_NUM];
  vals[V] = 1.;
  vals[X] = 0.;
  struct user_data udata = {10., 1.};
  EXIT_IF_0(AdamsSetYs0(adams_data, vals, EQUATIONS_NUM));
  EXIT_IF_0(AdamsSetX(adams_data, 0.));
  EXIT_IF_0(AdamsSetEquation(adams_data, RightSideV, V));
  EXIT_IF_0(AdamsSetEquation(adams_data, RightSideX, X));
  EXIT_IF_0(AdamsSetUserData(adams_data, &udata));
  EXIT_IF_0(AdamsSetStep(adams_data, STEP));
  EXIT_IF_0(AdamsCheck(adams_data));
  double t;
  do{
    AdamsStep(adams_data);
    t = AdamsGetX(adams_data);
    fprintf(a_res, "%.12g\t%.12g\t%.12g\n",
            t,
            AdamsGetY(adams_data, X),
            AdamsGetY(adams_data, V));
  }while(t <= 20.);
  AdamsFreeData(adams_data);
  fclose(a_res);
  return 1;
error:
  fclose(a_res);
  return 0;
}

int TestRK4(void){
  FILE * rkres = fopen("rk4.txt", "w");
  rk_data *rk4_data;
  EXIT_IF_0(RK4InitData(&rk4_data, EQUATIONS_NUM));
  double vals[EQUATIONS_NUM];
  vals[V] = 1.;
  vals[X] = 0.;
  struct user_data udata = {10., 1.};
  EXIT_IF_0(RK4SetYs0(rk4_data, vals, EQUATIONS_NUM));
  EXIT_IF_0(RK4SetX(rk4_data, 0.));
  EXIT_IF_0(RK4SetEquation(rk4_data, RightSideV, V));
  EXIT_IF_0(RK4SetEquation(rk4_data, RightSideX, X));
  EXIT_IF_0(RK4SetUserData(rk4_data, &udata));
  EXIT_IF_0(RK4SetStep(rk4_data, STEP));
  EXIT_IF_0(RK4Check(rk4_data));
  double t;
  do{
    RK4Step(rk4_data);
    t = RK4GetX(rk4_data);
    fprintf(rkres, "%.12g\t%.12g\t%.12g\n",
            t,
            RK4GetY(rk4_data, X),
            RK4GetY(rk4_data, V));
  }while(t <= 20.);
  RK4FreeData(rk4_data);
  fclose(rkres);
  return 1;
error:
  fclose(rkres);
  return 0;
}

int TestRK5(void){
  FILE *rkres = fopen("rk5.txt", "w");
  rk5_data *data;
  EXIT_IF_0(RK5InitData(&data, EQUATIONS_NUM));
  double vals[EQUATIONS_NUM];
  vals[V] = 1.;
  vals[X] = 0.;
  struct user_data udata = {10., 1.};
  EXIT_IF_0(RK5SetYs0(data, vals, EQUATIONS_NUM));
  EXIT_IF_0(RK5SetX(data, 0.));
  EXIT_IF_0(RK5SetEquation(data, RightSideV, V));
  EXIT_IF_0(RK5SetEquation(data, RightSideX, X));
  EXIT_IF_0(RK5SetUserData(data, &udata));
  EXIT_IF_0(RK5SetStep(data, STEP));
  EXIT_IF_0(RK5Check(data));
  double t;
  do{
    RK5Step(data);
    t = RK5GetX(data);
    fprintf(rkres, "%.12g\t%.12g\t%.12g\n",
            t,
            RK5GetY(data, X),
            RK5GetY(data, V));
  }while(t <= 20.);
  RK5FreeData(data);
  fclose(rkres);
  return 1;
error:
  fclose(rkres);
  return 0;
}

int TestAdams5(void){
  FILE *rkres = fopen("adams5.txt", "w");
  a5_data *data;
  EXIT_IF_0(Adams5InitData(&data, EQUATIONS_NUM));
  double vals[EQUATIONS_NUM];
  vals[V] = 1.;
  vals[X] = 0.;
  struct user_data udata = {10., 1.};
  EXIT_IF_0(Adams5SetYs0(data, vals, EQUATIONS_NUM));
  EXIT_IF_0(Adams5SetX(data, 0.));
  EXIT_IF_0(Adams5SetEquation(data, RightSideV, V));
  EXIT_IF_0(Adams5SetEquation(data, RightSideX, X));
  EXIT_IF_0(Adams5SetUserData(data, &udata));
  EXIT_IF_0(Adams5SetStep(data, STEP));
  EXIT_IF_0(Adams5Check(data));
  double t;
  do{
    Adams5Step(data);
    t = Adams5GetX(data);
    fprintf(rkres, "%.12g\t%.12g\t%.12g\n",
            t,
            Adams5GetY(data, X),
            Adams5GetY(data, V));
  }while(t <= 20.);
  Adams5FreeData(data);
  fclose(rkres);
  return 1;
error:
  fclose(rkres);
  return 0;
}

int main(int argc, char *argv[]){
  return !(TestAdams() && TestRK4() && TestRK5() && TestAdams5());
}
