#ifndef RK5_H
#define RK5_H

typedef double (*RK5RSFunc) (double const x,
                               double const *Y,
                               void *userdata);

typedef struct rk5_data_st rk5_data;

int RK5InitData(rk5_data **data, unsigned const eq_nums);
void RK5FreeData(rk5_data *data);
int RK5SetYs0(rk5_data *data, double const ys[],
                          unsigned const num);
int RK5SetY0(rk5_data *data, double const y, unsigned const index);
int RK5SetX(rk5_data *data, double const t);
int RK5SetStep(rk5_data *data, double const step);
int RK5SetEquation(rk5_data *data, RK5RSFunc func,
                     unsigned const index);
int RK5SetEquations(rk5_data *data, RK5RSFunc func[],
                      unsigned const num);
int RK5Check(rk5_data *data);
void RK5Step(rk5_data *data);
double RK5GetY(rk5_data *data, unsigned const num);
double *RK5GetYs(rk5_data *data);
double RK5GetX(rk5_data *data);
double RK5GetDY(rk5_data *data, unsigned const num);
int RK5SetUserData(rk5_data *data, void *userdata);

#endif //RK5_H
