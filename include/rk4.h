#ifndef RK4_H
#define RK4_H

typedef double (*RK4RSFunc) (double const x,
                               double const *Y,
                               void *userdata);

typedef struct rk4_data_st rk_data;

int RK4InitData(rk_data **data, unsigned const eq_nums);
void RK4FreeData(rk_data *data);
int RK4SetYs0(rk_data *data, double const ys[],
                          unsigned const num);
int RK4SetY0(rk_data *data, double const y, unsigned const index);
int RK4SetX(rk_data *data, double const t);
int RK4SetStep(rk_data *data, double const step);
int RK4SetEquation(rk_data *data, RK4RSFunc func,
                     unsigned const index);
int RK4SetEquations(rk_data *data, RK4RSFunc func[],
                      unsigned const num);
int RK4Check(rk_data *data);
void RK4Step(rk_data *data);
double RK4GetY(rk_data *data, unsigned const num);
double *RK4GetYs(rk_data *data);
double RK4GetX(rk_data *data);
double RK4GetDY(rk_data *data, unsigned const num);
int RK4SetUserData(rk_data *data, void *userdata);

#endif //RK4_H
