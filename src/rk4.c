#include <stdio.h>
#include "stdlib.h"
#include "rk4.h"

struct rk4_data_st{
  unsigned eq_num;
  double *y;
  double *f;
  double x;
  double h;
  RK4RSFunc *funcs;
  void *userdata;
};

#define EXIT_IF_NULL(POINTER) if( NULL == POINTER ){ goto error; }

int RK4InitData(rk_data **data, unsigned const eq_num){
  *data = calloc(1, sizeof(rk_data));
  EXIT_IF_NULL(*data);
  (*data)->eq_num = eq_num;
  (*data)->y = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->y);
  (*data)->funcs = calloc(eq_num, sizeof(RK4RSFunc));
  EXIT_IF_NULL((*data)->funcs);
  (*data)->f = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->f);
  return 1;
error:
  if(*data){
    free((*data)->funcs);
    free((*data)->y);
    free(*data);
    *data = NULL;
  }
  return 0;
}

void RK4FreeData(rk_data *data){
  if(data){
    free(data->f);
    free(data->funcs);
    free(data->y);
    free(data);
  }
}

int RK4SetYs0(rk_data *data, double const ys[],
                          unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->y[i] = ys[i];
  }
  return 1;
}

int RK4SetY0(rk_data *data, double const y, unsigned const index){
  if(!data || index >= data->eq_num){
    return 0;
  }
  data->y[index] = y;
  return 1;
}

int RK4SetX(rk_data *data, double const t){
  if(!data){
    return 0;
  }
  data->x = t;
  return 1;
}

int RK4SetStep(rk_data *data, double const step){
  if(!data){
    return 0;
  }
  data->h= step;
  return 1;
}

int RK4SetEquation(rk_data *data, RK4RSFunc func,
                     unsigned const index){
  if(!data || index >= data->eq_num){
      return 0;
    }
  data->funcs[index] = func;
  return 1;
}

int RK4SetEquations(rk_data *data, RK4RSFunc func[],
                      unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->funcs[i] = func[i];
  }
  return 1;
}

int RK4Check(rk_data *data){
  if(!data || !data->eq_num){
    fprintf(stderr, "%s\n", "RK4Check: Incorrect initialization.");
    return 0;
  }
  if(0. == data->h){
    fprintf(stderr, "%s\n", "RK4Check: Step must be greater then 0.");
    return 0;
  }
  for(unsigned i = 0; i< data->eq_num; i++){
    if(!data->funcs[i]){
      fprintf(stderr, "%s%d%s\n", "RK4Check: Right side functions for parameter number ", i, " not assigned.");
      return 0;
    }
  }
  return 1;
}

void RK4Step(rk_data *data){
  double k1[data->eq_num], k2[data->eq_num], k3[data->eq_num], k4[data->eq_num];
  double y[data->eq_num];
  double yn[data->eq_num];
  double h05 = data->h * 0.5;
  double t = data->x + h05;
  for(unsigned i = 0; i < data->eq_num; i++){
    k1[i] = data->funcs[i](data->x, data->y, data->userdata);
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
    data->f[i] = 1./6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    data->y[i] += data->h * data->f[i];
  }
}

double RK4GetY(rk_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->y[num];
}

double *RK4GetYs(rk_data *data){
  if(data){
    return data->y;
  }
  return NULL;
}

double RK4GetX(rk_data *data){
  if(data){
    return data->x;
  }
  return 0.;
}

double RK4GetDY(rk_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->f[num];
}

int RK4SetUserData(rk_data *data, void *userdata){
  if(data){
    data->userdata = userdata;
    return 1;
  }
  return 0;
}

