#include <stdio.h>
#include "stdlib.h"
#include "rk5.h"

struct rk5_data_st{
  unsigned eq_num;
  double *y;
  double *f;
  double x;
  double h;
  RK5RSFunc *funcs;
  void *userdata;
};

#define EXIT_IF_NULL(POINTER) if( NULL == POINTER ){ goto error; }

int RK5InitData(rk5_data **data, unsigned const eq_num){
  *data = calloc(1, sizeof(rk5_data));
  EXIT_IF_NULL(*data);
  (*data)->eq_num = eq_num;
  (*data)->y = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->y);
  (*data)->funcs = calloc(eq_num, sizeof(RK5RSFunc));
  EXIT_IF_NULL((*data)->funcs);
  (*data)->f = calloc(eq_num, sizeof(double));
  EXIT_IF_NULL((*data)->f);
  return 1;
error:
  if((*data)){
    free((*data)->funcs);
    free((*data)->y);
    free(*data);
    *data = NULL;
  }
  return 0;
}

void RK5FreeData(rk5_data *data){
  if(data){
    free(data->f);
    free(data->funcs);
    free(data->y);
    free(data);
  }
}

int RK5SetYs0(rk5_data *data, double const ys[],
                          unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->y[i] = ys[i];
  }
  return 1;
}

int RK5SetY0(rk5_data *data, double const y, unsigned const index){
  if(!data || index >= data->eq_num){
    return 0;
  }
  data->y[index] = y;
  return 1;
}

int RK5SetX(rk5_data *data, double const t){
  if(!data){
    return 0;
  }
  data->x = t;
  return 1;
}

int RK5SetStep(rk5_data *data, double const step){
  if(!data){
    return 0;
  }
  data->h= step;
  return 1;
}

int RK5SetEquation(rk5_data *data, RK5RSFunc func,
                     unsigned const index){
  if(!data || index >= data->eq_num){
      return 0;
    }
  data->funcs[index] = func;
  return 1;
}

int RK5SetEquations(rk5_data *data, RK5RSFunc func[],
                      unsigned const num){
  if(!data || num != data->eq_num){
    return 0;
  }
  for(unsigned i = 0; i<num; i++){
    data->funcs[i] = func[i];
  }
  return 1;
}

int RK5Check(rk5_data *data){
  if(!data || !data->eq_num){
    fprintf(stderr, "%s\n", "RK5Check: Incorrect initialization.");
    return 0;
  }
  if(0. == data->h){
    fprintf(stderr, "%s\n", "RK5Check: Step must be greater then 0.");
    return 0;
  }
  for(unsigned i = 0; i< data->eq_num; i++){
    if(!data->funcs[i]){
      fprintf(stderr, "%s%d%s\n", "RK5Check: Right side functions for parameter number ", i, " not assigned.");
      return 0;
    }
  }
  return 1;
}

static const double RK5_CONST[] = {1./24., 5./48., 27./56., 125./336.};

void RK5Step(rk5_data *data){
  double k1[data->eq_num], k2[data->eq_num], k3[data->eq_num],
         k4[data->eq_num], k5[data->eq_num], k6[data->eq_num];
  double y[data->eq_num];
  double yn[data->eq_num];
  for(unsigned i = 0; i < data->eq_num; i++){
    k1[i] = data->funcs[i](data->x, data->y, data->userdata);
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
    k6[i] = data->funcs[i](data->x + 1./5.*data->h,
                           y,
                           data->userdata);
    data->f[i] = RK5_CONST[0]*k1[i] + RK5_CONST[1]*k4[i] +
                 RK5_CONST[2]*k5[i] + RK5_CONST[3]*k6[i];
    data->y[i] += data->h*data->f[i];
  }
  data->x += data->h;
}

double RK5GetY(rk5_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->y[num];
}

double *RK5GetYs(rk5_data *data){
  if(data){
    return data->y;
  }
  return NULL;
}

double RK5GetX(rk5_data *data){
  if(data){
    return data->x;
  }
  return 0.;
}

double RK5GetDY(rk5_data *data, unsigned const num){
  if(!data || num >= data->eq_num){
    return 0.;
  }
  return data->f[num];
}

int RK5SetUserData(rk5_data *data, void *userdata){
  if(data){
    data->userdata = userdata;
    return 1;
  }
  return 0;
}


