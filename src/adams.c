#include <stdio.h>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "adams.h"

struct adams_data_st{
  unsigned eq_num;
  double *y;
  double *dy;
  double x;
  int boost_step;
  double h;
  AdamsRSFunc *funcs;
  double *f;
  void *userdata;
};

static const double kf[] = {55./24., -59./24., 37./24., -9./24.};
#define BOOST_STEPS 3

#define EXIT_IF_NULL(POINTER) if( NULL == POINTER ){ goto error; }

int AdamsInitData(a_data **data, unsigned const eq_num){
  *data = calloc(1, sizeof(a_data));
  EXIT_IF_NULL(*data);
  (*data)->eq_num = eq_num;
  (*data)->y = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->y);
  (*data)->funcs = calloc(eq_num, sizeof(AdamsRSFunc));
  EXIT_IF_NULL((*data)->funcs);
  (*data)->dy = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->dy);
  (*data)->f = calloc(eq_num*(BOOST_STEPS+1), sizeof(double));
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

void AdamsFreeData(a_data *data){
  if(data){
    free(data->f);
    free(data->funcs);
    free(data->y);
    free(data);
  }
}

int AdamsSetYs0(a_data *data, double const ys[],
                          unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->y[i] = ys[i];
  }
  return 1;
}

int AdamsSetY0(a_data *data, double const y, unsigned const index){
  if(!data || index >= data->eq_num){
    return 0;
  }
  data->y[index] = y;
  return 1;
}

int AdamsSetX(a_data *data, double const t){
  if(!data){
    return 0;
  }
  data->x = t;
  return 1;
}

int AdamsSetStep(a_data *data, double const step){
  if(!data){
    return 0;
  }
  data->h= step;
  data->boost_step = BOOST_STEPS;
  return 1;
}

int AdamsSetEquation(a_data *data, AdamsRSFunc func,
                     unsigned const index){
  if(!data || index >= data->eq_num){
      return 0;
    }
  data->funcs[index] = func;
  return 1;
}

int AdamsSetEquations(a_data *data, AdamsRSFunc func[],
                      unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->funcs[i] = func[i];
  }
  return 1;
}

int AdamsCheck(a_data *data){
  if(!data || !data->eq_num){
    fprintf(stderr, "%s\n", "AdamsCheck: Incorrect initialization.");
    return 0;
  }
  if(0. == data->h){
    fprintf(stderr, "%s\n", "AdamsCheck: Step must be greater then 0.");
    return 0;
  }
  for(unsigned i = 0; i< data->eq_num; i++){
    if(!data->funcs[i]){
      fprintf(stderr, "%s%d%s\n", "AdamsCheck: Right side functions for parameter number ", i, " not assigned.");
      return 0;
    }
  }
  return 1;
}

static void BoostRK4Step(a_data *data){
  double k1[data->eq_num], k2[data->eq_num], k3[data->eq_num], k4[data->eq_num];
  double y[data->eq_num];
  double yn[data->eq_num];
  double h05 = data->h * 0.5;
  double t = data->x + h05;
  (data->boost_step)--;
  for(unsigned i = 0; i < data->eq_num; i++){
    k1[i] = data->funcs[i](data->x, data->y, data->userdata);
    data->f[i*(BOOST_STEPS+1) + data->boost_step] = k1[i];
    y[i] = data->y[i] + h05*k1[i];
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k2[i] = data->funcs[i](t, y, data->userdata);
    yn[i] = data->y[i] + h05*k2[i];
  }
  for(unsigned i = 0; i < data->eq_num; i++){
    k3[i] = data->funcs[i](t, yn, data->userdata);
    y[i] = data->y[i] + data->h*k3[i];
  }
  data->x += data->h;
  for(unsigned i = 0; i < data->eq_num; i++){
    k4[i] = data->funcs[i](data->x, y, data->userdata);
    data->dy[i] = 1./6*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);
    data->y[i] += data->h*data->dy[i];
  }
}

static void MainAdamsStep(a_data *data){
  double y[data->eq_num];
  for(unsigned i = 0; i<data->eq_num; i++){
    unsigned j = i*(BOOST_STEPS + 1);
    memmove(&(data->f[j + 1]), &(data->f[j]), sizeof(double)*BOOST_STEPS);
    data->f[j] = data->funcs[i](data->x, data->y, data->userdata);
    data->dy[i] = kf[0]*data->f[j] +
                  kf[1]*data->f[j + 1] +
                  kf[2]*data->f[j + 2] +
                  kf[3]*data->f[j + 3];
    y[i] = data->y[i] + data->h*data->dy[i];
  }
  memcpy(data->y, y, sizeof(double)*data->eq_num);
  data->x += data->h;
}

void AdamsStep(a_data *data){
  if(data->boost_step){
    BoostRK4Step(data);
  } else {
    MainAdamsStep(data);
  }
}

double AdamsGetY(a_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->y[num];
}

double *AdamsGetYs(a_data *data){
  if(data){
    return data->y;
  }
  return NULL;
}

double AdamsGetX(a_data *data){
  if(data){
    return data->x;
  }
  return 0.;
}

double AdamsGetDY(a_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->dy[num];
}

int AdamsSetUserData(a_data *data, void *userdata){
  if(data){
    data->userdata = userdata;
    return 1;
  }
  return 0;
}
