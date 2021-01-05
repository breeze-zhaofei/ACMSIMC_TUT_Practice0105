#ifndef CONTROLLER_H
#define CONTROLLER_H

// #define VC_LOOP_CEILING 40

// struct PI_Reg{
//    double   Kp;
//    double   Ti;
//    double   Ki; // Ki = Kp/Ti*TS
//    double   i_state;
//    double   i_limit;
// };

// double PI(struct PI_Reg *r, double err);

struct ControllerForExperiment{

    double timebase;

    double ual;
    double ube;

    double rs;
    double rreq;
    double Lsigma;
    double alpha;
    double Lmu;
    double Lmu_inv;

    double Tload;
    double rpm_cmd;

    double Js;
    double Js_inv;

};
extern struct ControllerForExperiment CTRL;
void CTRL_init();


#endif