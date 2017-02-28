#include <stdio.h>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "adams5.h"

struct adams5_data_st{
  unsigned eq_num;
  double *y;
  double *dy;
  double x;
  int boost_step;
  double h;
  Adams5RSFunc *funcs;
  double *f;
  void *userdata;
};

static const double kf[] = {1901./720., -2774./720., 2616./720.,
                            -1274./720., 251./720.};

#define BOOST_STEPS 4

#define EXIT_IF_NULL(POINTER) if( NULL == POINTER ){ goto error; }

int Adams5InitData(a5_data **data, unsigned const eq_num){
  *data = calloc(1, sizeof(a5_data));
  EXIT_IF_NULL(*data);
  (*data)->eq_num = eq_num;
  (*data)->y = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->y);
  (*data)->funcs = calloc(eq_num, sizeof(Adams5RSFunc));
  EXIT_IF_NULL((*data)->funcs);
  (*data)->dy = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->dy);
  (*data)->f = calloc(eq_num*(BOOST_STEPS + 1), sizeof(double));
  EXIT_IF_NULL((*data)->f);
  return 1;
error:
  if(*data){
    free((*data)->dy);
    free((*data)->funcs);
    free((*data)->y);
    free(*data);
    *data = NULL;
  }
  return 0;
}

void Adams5FreeData(a5_data *data){
  if(data){
    free(data->f);
    free(data->dy);
    free(data->funcs);
    free(data->y);
    free(data);
  }
}

int Adams5SetYs0(a5_data *data, double const ys[],
                          unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->y[i] = ys[i];
  }
  return 1;
}

int Adams5SetY0(a5_data *data, double const y, unsigned const index){
  if(!data || index >= data->eq_num){
    return 0;
  }
  data->y[index] = y;
  return 1;
}

int Adams5SetX(a5_data *data, double const t){
  if(!data){
    return 0;
  }
  data->x = t;
  return 1;
}

int Adams5SetStep(a5_data *data, double const step){
  if(!data){
    return 0;
  }
  data->h= step;
  data->boost_step = BOOST_STEPS;
  return 1;
}

int Adams5SetEquation(a5_data *data, Adams5RSFunc func,
                     unsigned const index){
  if(!data || index >= data->eq_num){
      return 0;
    }
  data->funcs[index] = func;
  return 1;
}

int Adams5SetEquations(a5_data *data, Adams5RSFunc func[],
                      unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->funcs[i] = func[i];
  }
  return 1;
}

int Adams5Check(a5_data *data){
  if(!data || !data->eq_num){
    fprintf(stderr, "%s\n", "Adams5Check: Incorrect initialization.");
    return 0;
  }
  if(0. == data->h){
    fprintf(stderr, "%s\n", "Adams5Check: Step must be greater then 0.");
    return 0;
  }
  for(unsigned i = 0; i< data->eq_num; i++){
    if(!data->funcs[i]){
      fprintf(stderr, "%s%d%s\n", "Adams5Check: Right side functions for parameter number ", i, " not assigned.");
      return 0;
    }
  }
  return 1;
}

static const double RK5_CONST[] = {1./24., 5./48., 27./56., 125./336.};

static void BoostRK5Step(a5_data *data){
  double k1[data->eq_num], k2[data->eq_num], k3[data->eq_num],
         k4[data->eq_num], k5[data->eq_num], k6[data->eq_num];
  double y[data->eq_num];
  double yn[data->eq_num];
  (data->boost_step)--;
  for(unsigned i = 0; i < data->eq_num; i++){
    k1[i] = data->funcs[i](data->x, data->y, data->userdata);
    data->f[i*(BOOST_STEPS + 1) + data->boost_step] = k1[i];
    y[i] = data->y[i] + 0.5*data->h*k1[i];
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k2[i] = data->funcs[i](data->x + 0.5*data->h,
                           y,
                           data->userdata);
    yn[i] = data->y[i] + 0.25*data->h*(k1[i] + k2[i]);
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k3[i] = data->funcs[i](data->x + 0.5*data->h,
                           yn,
                           data->userdata);
    y[i] = data->y[i] + data->h*(2.*k3[i] - k2[i]);
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k4[i] = data->funcs[i](data->x + data->h,
                           y,
                           data->userdata);
    yn[i] = data->y[i] + 1./27.*data->h*(7.*k1[i] + 10.*k2[i]+k4[i]);
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k5[i] = data->funcs[i](data->x + 2./3.*data->h,
                           yn,
                           data->userdata);
    y[i] = data->y[i] + 1./625.*data->h*(28.*k1[i] - 125.*k2[i] +
                                         546.*k3[i] + 54.*k4[i] -
                                         378.*k5[i]);
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k6[i] = data->funcs[i](data->x + 1./5.*data->h, y, data->userdata);
    data->dy[i] = RK5_CONST[0]*k1[i] + RK5_CONST[1]*k4[i] +
                  RK5_CONST[2]*k5[i] + RK5_CONST[3]*k6[i];
    data->y[i] += data->h*data->dy[i];
  }
  data->x += data->h;
}

static void MainAdams5Step(a5_data *data){
  double y[data->eq_num];
  for(unsigned i = 0; i<data->eq_num; i++){
    unsigned j = i*(BOOST_STEPS + 1);
    memmove(&(data->f[j+1]), &(data->f[j]), sizeof(double)*BOOST_STEPS);
    data->f[j] = data->funcs[i](data->x, data->y, data->userdata);
    data->dy[i] = kf[0]*data->f[j] +
                  kf[1]*data->f[j + 1] +
                  kf[2]*data->f[j + 2] +
                  kf[3]*data->f[j + 3] +
                  kf[4]*data->f[j + 4];
    y[i] = data->y[i] + data->h*data->dy[i];
  }
  memcpy(data->y, y, sizeof(double)*data->eq_num);
  data->x += data->h;
}

void Adams5Step(a5_data *data){
  if(data->boost_step){
    BoostRK5Step(data);
  } else {
    MainAdams5Step(data);
  }
}

double Adams5GetY(a5_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->y[num];
}

double *Adams5GetYs(a5_data *data){
  if(data){
    return data->y;
  }
  return NULL;
}

double Adams5GetX(a5_data *data){
  if(data){
    return data->x;
  }
  return 0.;
}

double Adams5GetDY(a5_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->dy[num];
}

int Adams5SetUserData(a5_data *data, void *userdata){
  if(data){
    data->userdata = userdata;
    return 1;
  }
  return 0;
}
