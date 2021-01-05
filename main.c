#include "time.h"
#include "ACMSim.h"
#define NO_LOAD_TEST 1

struct InductionMachineSimulated IM;
void IM_init(){
    int i;
    for(i=0;i<5;++i){
        IM.x[i] = 0.0;
    }
    IM.rpm = 0.0;
    IM.rpm_cmd = 0.0;
    IM.rpm_deriv_cmd = 0.0;

    IM.iqs = 0.0;
    IM.ids = 0.0;

    IM.Tload = 0.0;
    IM.rpm_cmd = 0.0;
    IM.rpm_deriv_cmd = 0.0;

    // IM.Tem = 0.0;
    IM.Lmu    = 0.4482;
    IM.Lsigma = 0.0126;

    IM.rreq   = 1.69;
    IM.rs     = 3.04;

    IM.alpha  = IM.rreq / (IM.Lmu);
    IM.Lmu_inv= 1.0/IM.Lmu;

    IM.Js = 0.0636; // Awaya92 using im.omg
    IM.npp = 2;
    IM.mu_m = IM.npp/IM.Js;

    // IM.Ts  = MACHINE_TS;
    IM.Ts  = IM_TS;

    // IM.ial = 0.0;
    // IM.ibe = 0.0;
    IM.ual = 0.0;
    IM.ube = 0.0;
}


void rK5_dynamics(double t, double *x, double *fx){
        // electromagnetic model
        fx[2] = IM.rreq*x[0] - IM.alpha*x[2] - x[4]*x[3]; // flux-alpha
        fx[3] = IM.rreq*x[1] - IM.alpha*x[3] + x[4]*x[2]; // flux-beta
        fx[0] = (IM.ual - IM.rs*x[0] - fx[2])/IM.Lsigma; // current-alpha
        fx[1] = (IM.ube - IM.rs*x[1] - fx[3])/IM.Lsigma; // current-beta

        // mechanical model
        IM.Tem = IM.npp*(x[1]*x[2]-x[0]*x[3]);
        fx[4] = (IM.Tem -IM.Tload)*IM.mu_m; // elec. angular rotor speed
        // fx[5] = x[4];  

    // #if MACHINE_TYPE == INDUCTION_MACHINE
    //     // electromagnetic model
    //     fx[2] = ACM.rreq*x[0] - ACM.alpha*x[2] - x[4]*x[3]; // flux-alpha
    //     fx[3] = ACM.rreq*x[1] - ACM.alpha*x[3] + x[4]*x[2]; // flux-beta
    //     fx[0] = (ACM.ual - ACM.rs*x[0] - fx[2])/ACM.Lsigma; // current-alpha
    //     fx[1] = (ACM.ube - ACM.rs*x[1] - fx[3])/ACM.Lsigma; // current-beta

    //     // mechanical model
    //     ACM.Tem = ACM.npp*(x[1]*x[2]-x[0]*x[3]);
    //     fx[4] = (ACM.Tem - ACM.Tload)*ACM.mu_m; // elec. angular rotor speed
    //     fx[5] = x[4];                           // elec. angular rotor position

    // #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //     // electromagnetic model
    //     fx[0] = (ACM.ud - ACM.R * x[0] + x[2]*ACM.Lq*x[1]) / ACM.Ld; // current-d
    //     fx[1] = (ACM.uq - ACM.R * x[1] - x[2]*ACM.Ld*x[0] - x[2]*ACM.KE) / ACM.Lq; // current-q

    //     // mechanical model
    //     ACM.Tem = ACM.npp*(x[1]*ACM.KE + (ACM.Ld - ACM.Lq)*x[0]*x[1]);
    //     fx[2] = (ACM.Tem - ACM.Tload)*ACM.mu_m; // elec. angular rotor speed
    //     fx[3] = x[2];                           // elec. angular rotor position
    // #endif
}
void rK555_Lin(double t, double *x, double hs){
    // #if MACHINE_TYPE == INDUCTION_MACHINE
    //     #define NUMBER_OF_STATES 6
    // #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //     #define NUMBER_OF_STATES 4
    // #endif
    // #define NS NUMBER_OF_STATES

    // double k1[NS], k2[NS], k3[NS], k4[NS], xk[NS];
    // double fx[NS];
    double k1[5], k2[5], k3[5], k4[5], xk[5];
    double fx[5];
    int i,NS=5;

    rK5_dynamics(t, x, fx); // timer.t,
    for(i=0;i<NS;++i){        
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs, 
    for(i=0;i<NS;++i){        
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
    }
}
int machine_simulation(){
    rK555_Lin(CTRL.timebase, IM.x, IM.Ts);

        IM.ids  = IM.x[0];
        IM.iqs  = IM.x[1];
        IM.rpm  = IM.x[4] * 60 / (2 * M_PI * IM.npp);
        
        
        if(isNumber(IM.rpm))
            return false;
        else
            return true;        
        

    // API for explicit access
    // #if MACHINE_TYPE == INDUCTION_MACHINE
    //     ACM.ial    = ACM.x[0];
    //     ACM.ibe    = ACM.x[1];
    //     ACM.psi_al = ACM.x[2];
    //     ACM.psi_be = ACM.x[3];
    //     ACM.rpm    = ACM.x[4] * 60 / (2 * M_PI * ACM.npp);

    // #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //     ACM.theta_d = ACM.x[3];
    //     if(ACM.theta_d > M_PI){
    //         ACM.theta_d -= 2*M_PI;
    //     }else if(ACM.theta_d < -M_PI){
    //         ACM.theta_d += 2*M_PI; // 反转！
    //     }
    //     ACM.x[3] = ACM.theta_d;

    //     ACM.id  = ACM.x[0];
    //     ACM.iq  = ACM.x[1];
    //     ACM.ial = MT2A(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));
    //     ACM.ibe = MT2B(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));
    //     ACM.rpm = ACM.x[2] * 60 / (2 * M_PI * ACM.npp);
    // #endif

    // if(isNumber(ACM.rpm)){
    //     return false;
    // }else{
    //     printf("ACM.rpm is %g\n", ACM.rpm);
    //     return true;        
    // }
}
    /* Macro for External Access Interface */
// #define US(X) im.us[X]
// #define IS(X) im.is[X]
// #define US_C(X) im.us_curr[X]
// #define IS_C(X) im.is_curr[X]
// #define US_P(X) im.us_prev[X]
// #define IS_P(X) im.is_prev[X]
void measurement(){
    US_C(0) = CTRL.ual;
    US_C(1) = CTRL.ube;
    US_P(0) = US_C(0);
    US_P(1) = US_C(1);

    IS_C(0) = IM.ids;
    IS_C(1) = IM.iqs;

    im.omg = IM.x[4];

    // #if MACHINE_TYPE == INDUCTION_MACHINE
    //     IS_C(0) = ACM.ial;
    //     IS_C(1) = ACM.ibe;
    //     im.omg = ACM.x[4];
    //     im.theta_r = ACM.x[5];
    // #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //     IS_C(0) = ACM.ial;
    //     IS_C(1) = ACM.ibe;
    //     sm.omg = ACM.x[2];
    //     sm.theta_d = ACM.x[3];
    //     sm.theta_r = sm.theta_d;
    // #endif
}
void inverter_model(){

    IM.ual = CTRL.ual;
    IM.ube = CTRL.ube;

    // #if MACHINE_TYPE == INDUCTION_MACHINE
    //     ACM.ual = CTRL.ual;
    //     ACM.ube = CTRL.ube;
    // #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //     ACM.ual = CTRL.ual;
    //     ACM.ube = CTRL.ube;
    //     ACM.ud = AB2M(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
    //     ACM.uq = AB2T(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
    // #endif
}

int main(){
    printf("NUMBER_OF_LINES: %d\n\n", NUMBER_OF_LINES);

    /* Initialization */
    // Machine_init();
    IM_init();
    CTRL_init();
    // acm_init();
    // ob_init();

    FILE *fw;
    fw = fopen("algorithm.dat", "w");
    // write_header_to_file(fw);

    /* MAIN LOOP */
    clock_t  begin, end;
    begin = clock();
    int _; // _ for the outer iteration
    int dfe=0; // dfe for down frequency execution
    for(_=0;_<NUMBER_OF_LINES;++_){

        /* Command and Load Torque */
        IM.rpm_cmd = 50;
        IM.Tload = 0;
        // ACM.Tload = 0;
        // cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 1500); // timebase, instant, interval, rpm_cmd
        // if(CTRL.timebase>10){
        //     ACM.rpm_cmd = -250;
        // }else if(CTRL.timebase>5){
        //     ACM.Tload = 5;
        // }else{
        //     ACM.rpm_cmd = -50;
        //     ACM.Tload = 1;
        // }

        /* Simulated ACM */
        if(machine_simulation()){ 
            printf("Break the loop.\n");
            break;
        }

        if(++dfe==DOWN_FREQ_EXE){
            dfe = 0;

            /* Time */
            CTRL.timebase += TS;

            measurement();

            // observation();

            write_data_to_file(fw);

        #if NO_LOAD_TEST == TRUE
            #define VF_RATIO 18 //18.0 // 8 ~ 18 shows saturated phenomenon
            double freq = 2; // 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
            double volt = VF_RATIO*freq;
            CTRL.ual = volt*cos(2*M_PI*freq*CTRL.timebase);
            CTRL.ube = volt*sin(2*M_PI*freq*CTRL.timebase);
        #else
            // control(ACM.rpm_cmd, 0);
            control();
        #endif
         
        }

        inverter_model();
    }
    end = clock(); printf("The simulation in C costs %g sec.\n", (double)(end - begin)/CLOCKS_PER_SEC);
    fclose(fw);

    /* Fade out */
    system("python ./ACMPlot.py"); 
    // getch();
    // system("pause");
    // system("exit");
    return 0; 
}


void write_data_to_file(FILE *fw){
    // static int bool_animate_on = false;
    static int j=0,jj=0; // j,jj for down sampling

    // if(CTRL.timebase>20)
    if(++j == 10)
    {
        j=0;
        // 数目必须对上，否则ACMAnimate会失效，但是不会影响ACMPlot
        fprintf(fw, "%g,%g,%g,%g,%g\n",
        IM.x[0], IM.x[1], IM.x[2], IM.x[3], IM.x[4] );

    }


    // {
    //     if(++j == DOWN_SAMPLE)
    //     {
    //         j=0;
    //         #if MACHINE_TYPE == INDUCTION_MACHINE
    //             // 数目必须对上，否则ACMAnimate会失效，但是不会影响ACMPlot
    //             fprintf(fw, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
    //                     ACM.x[0], ACM.x[1], ACM.x[2], ACM.x[3], ACM.x[4], ACM.Tem, 
    //                     CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.iMs_cmd, CTRL.iMs, CTRL.iTs_cmd, CTRL.iTs,
    //                     ob.psi_mu_al, ob.tajima.omg*RAD_PER_SEC_2_RPM, ACM.rpm
    //                     );
    //         #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    //             fprintf(fw, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
    //                     ACM.x[0], ACM.x[1], ACM.x[2], ACM.x[3],
    //                     CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.iMs_cmd, CTRL.iMs, CTRL.iTs_cmd, CTRL.iTs
    //                     );
    //         #endif
    //     }
    // }

    // if(bool_animate_on==false){
    //     bool_animate_on = true;
    //     printf("Start ACMAnimate\n");
    //     system("start python ./ACMAnimate.py"); 
    // }
}

int isNumber(double x){
    // This looks like it should always be true, 
    // but it's false if x is a NaN (1.#QNAN0).
    return (x == x); 
    // see https://www.johndcook.com/blog/IEEE_exceptions_in_cpp/ cb: https://stackoverflow.com/questions/347920/what-do-1-inf00-1-ind00-and-1-ind-mean
}




