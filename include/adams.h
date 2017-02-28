#ifndef ADAMS_H
#define ADAMS_H

typedef double (*AdamsRSFunc) (double const x,
                               double const *Y,
                               void *userdata);

typedef struct adams_data_st a_data;

int AdamsInitData(a_data **data, unsigned const eq_nums);
void AdamsFreeData(a_data *data);
int AdamsSetYs0(a_data *data, double const ys[],
                          unsigned const num);
int AdamsSetY0(a_data *data, double const y, unsigned const index);
int AdamsSetX(a_data *data, double const t);
int AdamsSetStep(a_data *data, double const step);
int AdamsSetEquation(a_data *data, AdamsRSFunc func,
                     unsigned const index);
int AdamsSetEquations(a_data *data, AdamsRSFunc func[],
                      unsigned const num);
int AdamsCheck(a_data *data);
void AdamsStep(a_data *data);
double AdamsGetY(a_data *data, unsigned const num);
double *AdamsGetYs(a_data *data);
double AdamsGetX(a_data *data);
double AdamsGetDY(a_data *data, unsigned const num);
int AdamsSetUserData(a_data *data, void *userdata);

#endif //ADAMS_H
