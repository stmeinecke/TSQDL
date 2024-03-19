class DDEintegrator{
  //this class provides Runge-Kutta methods to solve delay differential equations. As of now, there are fixed stepsize first order (Euler), second order and fourth order methods which all take the same arguments:
  //1. A std::function "derivs" that defines the dynamical system "derivs". Can be a Lambda-function defined inside the the main function. The derivs function takes the arguments: 1. current time (double), 2. current state (pointer to a vars object), 3. history vector (pointer to a vars_vec object), 4. derivatives which are to be saved in a (pointer to a vars object), 5. system parameters (pointer to a parameters object).
  //2. A std::function "algs" that performs algebraic equation on the dynamical variables. Can be a Lambda-function defined inside the the main function. Takes the arguments: current time, current system state, system history and system parameters. Is carried out after the system has been advanced after a step. Caution! Therefore, the system is at a different state at_t0 compared to th RK calculations. Use the provided function algs_empty if no algebraic equations exist.
  //3. A std::function "noise" that calculates a noise parameter for each RK step
  //4. A pointer to a vars_vec history vector "X" of the dynamical system, which also stores the current state of the system
  //5. A pointer to a parameters object "paras" that contains all the relevant system parameters
  //6. A pointer to the dynamical time "t", which will be propagated during the integration
  //7. The integration time "T". Not the the integration is performed from t to t+T and not from t to T
  //8. The timestep "dt"
  //9. A std::function "outputfunc" that is executed before every time step to allow the output of the current state. Can be a Lambda-function defined inside the the main function

  public:
    void DDE_euler(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    void DDE_euler(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    void DDE_RK2(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    void DDE_RK2(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    void DDE_RK4(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
    void DDE_RK4(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc);
  private:
};

void DDEintegrator::DDE_euler(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for euler method
  vars dXdt;
  dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 0;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate derivatives at the current position
    derivs(Xhist->at_t0_rPtr(),Xhist, &dXdt, paras);
    dXdt.mult(&dXdt, dt);
    
    //calculate next step
    Xhist->at_t0_rPtr()->add(&dXdt,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
        
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
  }
}

void DDEintegrator::DDE_RK2(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK2 method
  vars dXdt, k1, k2;
  k1.setTo(0.0), k2.setTo(0.0), dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 1;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate trial step across dt and save it to k1
    derivs(Xhist->at_t0_rPtr(), Xhist, &dXdt, paras);
    dXdt.mult(&k1, dt/2.0);
    Xhist->at_t0_rPtr()->add(&k1, &k1);
    
    //calculate the step k2 using the trial step k1
    Xhist->tau_offset = 0.5*dt; //adjust tau_offset for the intermediate step
    algs(&k1, Xhist, paras);
    derivs(&k1, Xhist, &dXdt, paras);
    Xhist->tau_offset = 0.0; //reset tau_offset
    dXdt.mult(&k2,dt);
    
    //calculate the next state using k2 and save it to the history vector
    Xhist->at_t0_rPtr()->add(&k2,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
        
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
  }
}

void DDEintegrator::DDE_RK4(std::function<void (vars*, vars_vec*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK4 method
  vars k1, k2, k3, k4, xtmp;
  k1.setTo(0.0), k2.setTo(0.0), k3.setTo(0.0), k4.setTo(0.0), xtmp.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 3;

  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate next step -> results are written to Xhist->at_next = one step in the future
    //calculate RK4 indvividual components ki that add to x as in x_{n+1} = x_{n} + dt*(k'1/6 + k'2/3 + k3/3 + k4/6) where
    
    //first step
    derivs(Xhist->at_t0_rPtr(),Xhist,&k1,paras);
    k1.mult(&k1,dt/2.0); //k1 = dt/2*f(x,t,paras) = k'1*0.5
    
    //second step
    Xhist->at_t0_rPtr()->add(&k1,&xtmp); //xtmp = x + k1
    Xhist->tau_offset = 0.5*dt;
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k2,paras); 
    k2.mult(&k2,dt/2.0); //k2 = dt/2*f(xtmp=x + k1,t + dt/2,paras) = k'2*0.5
    
    //third step
    Xhist->at_t0_rPtr()->add(&k2,&xtmp); //xtmp = x + k2
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k3,paras);
    k3.mult(&k3,dt); //k3 = dt*f(xtemp = x + k2, t + dt/2,paras)
    
    //fourth step
    Xhist->at_t0_rPtr()->add(&k3,&xtmp); //xtemp = x + k3
    Xhist->tau_offset = dt;
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k4,paras);
    k4.mult(&k4,dt); //k4 = dt*f(xtemp = x + k3, t + dt,paras)
    //advance x by  k1/3 + k2/1.5 + k3/3 + k4/6
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
	  ((double*)Xhist->at_nxt_rPtr())[i] = ((double*)Xhist->at_t0_rPtr())[i]
					+ ((double*)&k1)[i] / 3.0
					+ ((double*)&k2)[i] / 1.5
					+ ((double*)&k3)[i] / 3.0
					+ ((double*)&k4)[i] / 6.0;
    }
    
    Xhist->tau_offset = 0.0; //reset tau_offset
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
    
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
    
  }
}

/////////////////////
//with vars_vec_wdX
/////////////////////

void DDEintegrator::DDE_euler(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for euler method
  vars dXdt;
  dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 0;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate derivatives at the current position
    derivs(Xhist->at_t0_rPtr(),Xhist, &dXdt, paras);
    dXdt.mult(&dXdt, dt);
    
    //calculate next step
    Xhist->at_t0_rPtr()->add(&dXdt,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
        
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
  }
}

void DDEintegrator::DDE_RK2(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK2 method
  vars dXdt, k1, k2;
  k1.setTo(0.0), k2.setTo(0.0), dXdt.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 1;
  
  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate trial step across dt and save it to k1
    derivs(Xhist->at_t0_rPtr(), Xhist, &dXdt, paras);
    dXdt.mult(&k1, dt/2.0);
    Xhist->at_t0_rPtr()->add(&k1, &k1);
    
    //calculate the step k2 using the trial step k1
    Xhist->tau_offset = 0.5*dt; //adjust tau_offset for the intermediate step
    algs(&k1, Xhist, paras);
    derivs(&k1, Xhist, &dXdt, paras);
    Xhist->tau_offset = 0.0; //reset tau_offset
    dXdt.mult(&k2,dt);
    
    //calculate the next state using k2 and save it to the history vector
    Xhist->at_t0_rPtr()->add(&k2,Xhist->at_nxt_rPtr());
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn+=1;
        
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
  }
}


void DDEintegrator::DDE_RK4(std::function<void (vars*, vars_vec_wdX*, vars*, parameters*)> derivs, std::function<void (vars*, vars_vec_wdX*, parameters*)> algs, std::function<void (parameters*)> noise, vars_vec_wdX* Xhist, parameters* paras, unsigned long long int* tn, const double T, const double dt, std::function<void (vars*, unsigned long long int, unsigned long long int)> outputfunc){
  //initialize vars for RK4 method
  vars k1, k2, k3, k4, xtmp;
  k1.setTo(0.0), k2.setTo(0.0), k3.setTo(0.0), k4.setTo(0.0), xtmp.setTo(0.0);
  
  //set history vector interpolation order
  Xhist->IntPol_order = 3;

  //calc upper integration time limit
  unsigned long long int tn_final = *tn + (unsigned long long int)(T/dt);

  //integration loop
  while(*tn <= tn_final){    
    //output to vector/file...
    outputfunc(Xhist->at_t0_rPtr(), *tn, tn_final);
    
    //calc noise
    noise(paras);
    
    //calculate next step -> results are written to Xhist->at_next = one step in the future
    //calculate RK4 indvividual components ki that add to x as in x_{n+1} = x_{n} + dt*(k'1/6 + k'2/3 + k3/3 + k4/6) where
    
    //first step
    derivs(Xhist->at_t0_rPtr(),Xhist,&k1,paras);
    //write derivative to derivative history vector
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
      ((double*)Xhist->dX_at_t0_rPtr())[i] = ((double*)&k1)[i];
    }
    //mult k1 with dt/2
    k1.mult(&k1,dt/2.0); //k1 = dt/2*f(x,t,paras) = k'1*0.5
    
    //second step
    Xhist->at_t0_rPtr()->add(&k1,&xtmp); //xtmp = x + k1
    Xhist->tau_offset = 0.5*dt;
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k2,paras); 
    k2.mult(&k2,dt/2.0); //k2 = dt/2*f(xtmp=x + k1,t + dt/2,paras) = k'2*0.5
    
    //third step
    Xhist->at_t0_rPtr()->add(&k2,&xtmp); //xtmp = x + k2
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k3,paras);
    k3.mult(&k3,dt); //k3 = dt*f(xtemp = x + k2, t + dt/2,paras)
    
    //fourth step
    Xhist->at_t0_rPtr()->add(&k3,&xtmp); //xtemp = x + k3
    Xhist->tau_offset = dt;
    algs(&xtmp, Xhist, paras);
    derivs(&xtmp,Xhist,&k4,paras);
    k4.mult(&k4,dt); //k4 = dt*f(xtemp = x + k3, t + dt,paras)
    //advance x by  k1/3 + k2/1.5 + k3/3 + k4/6
    for(std::size_t i=0; i<sizeof(vars)/sizeof(double); i+=1 ){
      ((double*)Xhist->at_nxt_rPtr())[i] = ((double*)Xhist->at_t0_rPtr())[i]
                                        + ((double*)&k1)[i] / 3.0
                                        + ((double*)&k2)[i] / 1.5
                                        + ((double*)&k3)[i] / 3.0
                                        + ((double*)&k4)[i] / 6.0;
    }
    
    Xhist->tau_offset = 0.0; //reset tau_offset
    
    //update t0 in the history vector and increment t by dt
    Xhist->incr_t0();
    *tn += 1;
    
    //after step algebraic functions -> RK calculations have been already carried out and t0 has been advanced -> Xhist->at_t0 is one step further compared to the RK calculations
    algs(Xhist->at_t0_rPtr(), Xhist, paras);
    
  }
}


auto outputfunc_empty = [](vars* X, unsigned long long int tn, unsigned long long int tn_final){};
auto algs_empty = [](vars* X, vars_vec* Xhist, parameters* p){};
auto noise_empty = [](parameters* p){};




