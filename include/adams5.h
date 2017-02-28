#ifndef ADAMS5_H
#define ADAMS5_H

typedef double (*Adams5RSFunc) (double const x,
                               double const *Y,
                               void *userdata);

typedef struct adams5_data_st a5_data;

int Adams5InitData(a5_data **data, unsigned const eq_nums);
void Adams5FreeData(a5_data *data);
int Adams5SetYs0(a5_data *data, double const ys[],
                          unsigned const num);
int Adams5SetY0(a5_data *data, double const y, unsigned const index);
int Adams5SetX(a5_data *data, double const t);
int Adams5SetStep(a5_data *data, double const step);
int Adams5SetEquation(a5_data *data, Adams5RSFunc func,
                     unsigned const index);
int Adams5SetEquations(a5_data *data, Adams5RSFunc func[],
                      unsigned const num);
int Adams5Check(a5_data *data);
void Adams5Step(a5_data *data);
double Adams5GetY(a5_data *data, unsigned const num);
double *Adams5GetYs(a5_data *data);
double Adams5GetX(a5_data *data);
double Adams5GetDY(a5_data *data, unsigned const num);
int Adams5SetUserData(a5_data *data, void *userdata);

#endif //ADAMS5_H
