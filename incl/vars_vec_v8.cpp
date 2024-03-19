//what?:
//vars_vec is a wrapper class for a vector of vars to be used as a history function for the DDE solver. It keeps track of the most recent entry t0 and provides methods a access or interpolate variables at a given delay tau < tau_max.

//the methods at_R (for real valued variables) and at_C (for complex valued variables) select an interpolation order based upon the "IntPol_order" value, which is set by the DDEsolver. Thus, they are recommended to be used in the formulation of the DDEs.


class vars_vec{
  //still a bit messy around here, but it works
  
  private:

  public:
    double dt, tau_max;
    unsigned long long int tn;
    unsigned int hist_length;
    unsigned int t0;
    bool dt_divides_tau = false;
    
    std::vector<vars> Xt;
    double tau_offset; //used to adjust tau for during in between steps of Runge-Kutta methods
    int IntPol_order; //set by the solver to select the appropriate delay vector interpolation
  
  vars_vec(double tau_max, double dt){
    //constructor creates a vector of appropriate length such that hist_length * dt >= tau_max
    this->tau_max = tau_max;
    this->dt = dt;
    this->hist_length = ceil(tau_max/dt) + 2; //ceil makes sure the vector is big enough, +2 in order to keep one additional point for the 3rd order interpolation and another to also save the current state
    Xt.resize(hist_length);
    this->t0 = hist_length-1; //current state is stored in the last element of the vector
    this->tau_offset = 0.0;
    this->IntPol_order = 1;
    this->tn = 0;
  }
  
  void setToCnst(double s){
    //sets all variables of all vars to the input s
    for (auto & v : this->Xt) {
      v.setTo(s);
    }
  }
  
  void setToCnst(double s, int index){
    //sets variable v of all vars to the input s
    for (auto & v : this->Xt) {
      ((double*)&(v))[index] = s;
    }
  }
  
  void save(std::string filename){
    std::ofstream binout_Xhist (filename, std::ios::out | std::ios::binary);
    
    for(std::size_t i=0; i<sizeof(int); i++) binout_Xhist << ((unsigned char*)(&(hist_length)))[i];
    for(std::size_t i=0; i<sizeof(int); i++) binout_Xhist << ((unsigned char*)(&(t0)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(dt)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(tau_max)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(tn)))[i];
    
    for(unsigned int k = 0; k < hist_length; k++){
      for(unsigned int l = 0; l < sizeof(vars); l++){
        binout_Xhist << ((unsigned char*)(&this->Xt[k]))[l];
      }
    }
    
//     std::cout << "History saved" << std::endl;
    
    binout_Xhist.close();
  }
  
  void load(std::string filename){
    std::ifstream binin_Xhist (filename, std::ios::in | std::ios::binary);
    if(binin_Xhist.good()){
    
      for(std::size_t i=0; i<sizeof(int); i++) ((unsigned char*)(&(hist_length)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(int); i++) ((unsigned char*)(&(t0)))[i] = binin_Xhist.get();
      this->Xt.resize(this->hist_length);
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(dt)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(tau_max)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(tn)))[i] = binin_Xhist.get();
      
      for(unsigned int k = 0; k < hist_length; k++){
        for(unsigned int l = 0; l < sizeof(vars); l++){
          ((unsigned char*)(&this->Xt[k]))[l] = binin_Xhist.get();
        }
      }
      
//       std::cout << "History loaded" << std::endl;
    }
    else std::cout << "Error: History file not found" << std::endl;
    
    binin_Xhist.close();
  }
  
  void incr_t0(){
    //increment current time index t0 by one
    this-> t0 = (this->t0+1)%this->hist_length;
  }

  vars * at_t0_rPtr(){
    //returns pointer to vars at current time t0
    return &(this->Xt[this->t0]);
  }
  
  vars * at_nxt_rPtr(){
    //returns pointer to the newest (next) vars, which previously was the oldest vars in the history
    return &(this->Xt[(this->t0 + 1)%this->hist_length]);
  }

  vars * at_Rnd_rPtr(double tau){
    //return pointer to vars at time t0-tau where tau/dt is rounded to the closest integer
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return nullptr;
    }
    else{
      return &(this->Xt[(this->t0 + this->hist_length - (int)round((tau-this->tau_offset)/this->dt)) % this->hist_length]);
    }
  }
  
  //////////////////////////////////////
  // save and general version
  //////////////////////////////////////
  
  double at_R(double tau, int index){
    if(this->IntPol_order == 3){
      if(this->dt_divides_tau == true){
        if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
          return at_IntPol_3rd_R_hlfstp(tau, index);
        }
        else{
          return ((double*)this->at_Rnd_rPtr(tau))[index];
        }
      }
      else return at_IntPol_3rd_R(tau,index);
    }
    else if(this->IntPol_order == 1) return at_IntPol_Lin_R(tau,index);
    else if(this->IntPol_order == 0) return ((double*)this->at_Rnd_rPtr(tau))[index];
    else{
      std::cout << "Error: interpolation order does not exist" << std::endl;
      return 0.0;
    }
  }
  
  std::complex<double> at_C(double tau, int index){
    if(this->IntPol_order == 3){
      if(this->dt_divides_tau == true){
        if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
          return at_IntPol_3rd_C_hlfstp(tau, index);
        }
        else{
          return std::complex<double>(((double*)this->at_Rnd_rPtr(tau))[index],((double*)this->at_Rnd_rPtr(tau))[index+1]);
        }
      }
      else return at_IntPol_3rd_C(tau,index);
    }
    else if(this->IntPol_order == 1) return at_IntPol_Lin_C(tau,index);
    else if(this->IntPol_order == 0) return std::complex<double>(((double*)this->at_Rnd_rPtr(tau))[index],((double*)this->at_Rnd_rPtr(tau))[index+1]);
    else{
      std::cout << "Error: interpolation order does not exist" << std::endl;
      return std::complex<double> (0.0,0.0);
    }
  }
  
  //////////////////////////////////////
  // fast 3rd interpol with halfsteps
  //////////////////////////////////////
  
  double at_R_fast3rd(double tau, int index){
    if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
      return at_IntPol_3rd_R_hlfstp(tau, index);
    }
    else{
      return ((double*)this->at_Rnd_rPtr(tau))[index];
    }
  }
  
  std::complex<double> at_C_fast3rd(double tau, int index){
    if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
      return at_IntPol_3rd_C_hlfstp(tau, index);
    }
    else{
      int rind = (this->t0 + this->hist_length - (int)round((tau-this->tau_offset)/this->dt)) % this->hist_length;
      return std::complex<double>(((double*)(&(this->Xt[rind])))[index],((double*)(&(this->Xt[rind])))[index+1]);
    }
  }
  
  //////////////////////////////////////
  // direct 1st order interpol
  //////////////////////////////////////
  
  double at_IntPol_Lin_R(double tau, int index){
    //returns a double of the indexed variable of vars which is linearly interpolated at the time t0-tau
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return 0.0;
    }
    else{
      int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
      int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
      double rel_dt = (tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt);
      return ((double*)&(this->Xt.at(lower_hist)))[index] * (1-rel_dt) + ((double*)&(this->Xt.at(upper_hist)))[index] * rel_dt;
    }
  }
  
  std::complex<double> at_IntPol_Lin_C(double tau, int index){
    //returns a complex of the indexed complex variable of vars which is linearly interpolated at the time t0-tau
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return std::complex<double> (0.0,0.0);
    }
    else{
      int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
      int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
      double rel_dt = (tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt);
      return std::complex<double>(((double*)&(this->Xt.at(lower_hist)))[index] * (1-rel_dt) + ((double*)&(this->Xt.at(upper_hist)))[index] * rel_dt, ((double*)&(this->Xt.at(lower_hist)))[index+1] * (1-rel_dt) + ((double*)&(this->Xt.at(upper_hist)))[index+1] * rel_dt);
    }
  }
  
  //////////////////////////////////////
  // direct 3rd order interpol
  //////////////////////////////////////
  
  double at_IntPol_3rd_R(double tau, int index){
    //returns a double of the indexed variable of vars which is interpolated to third order using Neville's algorithm at the time t0-tau
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return 0;
    }
    else{
      int lower_hist = (this->t0 + 2*this->hist_length - (int)floor((tau-this->tau_offset)/this->dt) - 2) % this->hist_length;
      double rel_tau = 2*this->dt - this->dt*((tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt));
      double p00 = ((double*)&(this->Xt.at((lower_hist+0)%this->hist_length)))[index];
      double p11 = ((double*)&(this->Xt.at((lower_hist+1)%this->hist_length)))[index];
      double p22 = ((double*)&(this->Xt.at((lower_hist+2)%this->hist_length)))[index];
      double p33 = ((double*)&(this->Xt.at((lower_hist+3)%this->hist_length)))[index];
      
      double p01 = ((1.0*dt-rel_tau)*p00 + (rel_tau-0.0*dt)*p11)/this->dt;
      double p12 = ((2.0*dt-rel_tau)*p11 + (rel_tau-1.0*dt)*p22)/this->dt;
      double p23 = ((3.0*dt-rel_tau)*p22 + (rel_tau-2.0*dt)*p33)/this->dt;
      
      double p02 = ((2.0*dt-rel_tau)*p01 + (rel_tau-0.0*dt)*p12)/(2.0*this->dt);
      double p13 = ((3.0*dt-rel_tau)*p12 + (rel_tau-1.0*dt)*p23)/(2.0*this->dt);
      
      double p03 = ((3.0*dt-rel_tau)*p02 + (rel_tau-0.0*dt)*p13)/(3.0*this->dt);
      
      return p03;
    }
  }
  
  std::complex<double> at_IntPol_3rd_C(double tau, int index){
    //returns a double of the indexed variable of vars which is interpolated to third order using Neville's algorithm at the time t0-tau
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return 0;
    }
    else{
      
      int lower_hist = (this->t0 + 2*this->hist_length - (int)floor((tau-this->tau_offset)/this->dt) - 2) % this->hist_length;
      double rel_tau = 2*this->dt - this->dt*((tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt));
      //interpolate real part
      double p00 = ((double*)&(this->Xt.at((lower_hist+0)%this->hist_length)))[index];
      double p11 = ((double*)&(this->Xt.at((lower_hist+1)%this->hist_length)))[index];
      double p22 = ((double*)&(this->Xt.at((lower_hist+2)%this->hist_length)))[index];
      double p33 = ((double*)&(this->Xt.at((lower_hist+3)%this->hist_length)))[index];
      
      double p01 = ((1.0*dt-rel_tau)*p00 + (rel_tau-0.0*dt)*p11)/this->dt;
      double p12 = ((2.0*dt-rel_tau)*p11 + (rel_tau-1.0*dt)*p22)/this->dt;
      double p23 = ((3.0*dt-rel_tau)*p22 + (rel_tau-2.0*dt)*p33)/this->dt;
      
      double p02 = ((2.0*dt-rel_tau)*p01 + (rel_tau-0.0*dt)*p12)/(2.0*this->dt);
      double p13 = ((3.0*dt-rel_tau)*p12 + (rel_tau-1.0*dt)*p23)/(2.0*this->dt);
      
      double p03R = ((3.0*dt-rel_tau)*p02 + (rel_tau-0.0*dt)*p13)/(3.0*this->dt);
      
      //interpolate imaginary part
      p00 = ((double*)&(this->Xt.at((lower_hist+0)%this->hist_length)))[index+1];
      p11 = ((double*)&(this->Xt.at((lower_hist+1)%this->hist_length)))[index+1];
      p22 = ((double*)&(this->Xt.at((lower_hist+2)%this->hist_length)))[index+1];
      p33 = ((double*)&(this->Xt.at((lower_hist+3)%this->hist_length)))[index+1];
      
      p01 = ((1.0*dt-rel_tau)*p00 + (rel_tau-0.0*dt)*p11)/this->dt;
      p12 = ((2.0*dt-rel_tau)*p11 + (rel_tau-1.0*dt)*p22)/this->dt;
      p23 = ((3.0*dt-rel_tau)*p22 + (rel_tau-2.0*dt)*p33)/this->dt;
      
      p02 = ((2.0*dt-rel_tau)*p01 + (rel_tau-0.0*dt)*p12)/(2.0*this->dt);
      p13 = ((3.0*dt-rel_tau)*p12 + (rel_tau-1.0*dt)*p23)/(2.0*this->dt);
      
      double p03C = ((3.0*dt-rel_tau)*p02 + (rel_tau-0.0*dt)*p13)/(3.0*this->dt);
      
      return std::complex<double>(p03R,p03C);
    }
  }
  
  //////////////////////////////////////
  // direct third order with first order for small tau
  //////////////////////////////////////
  
  double at_R_3rd_wSmallTau(double tau, int index){
    if(tau < 2.0*dt){
      return at_IntPol_Lin_R(tau,index);
    }
    else{
      return at_IntPol_3rd_R(tau,index);
    }
  }
  
  std::complex<double> at_C_3rd_wSmallTau(double tau, int index){
    if(tau < 2.0*dt){
      return at_IntPol_Lin_C(tau,index);
    }
    else{
      return at_IntPol_3rd_C(tau,index);
    }
  }
  
  //////////////////////////////////////
  // only to be used to calc halfsteps
  //////////////////////////////////////
  
  double at_IntPol_3rd_R_hlfstp(double tau, int index){
  //returns a double of the indexed variable of vars which is interpolated to third order using Neville's algorithm at the time t0-tau
    int lower_hist = (this->t0 + 2*this->hist_length - (int)floor((tau-this->tau_offset)/this->dt) - 2) % this->hist_length;
    //interpolate according to p03 = (-p00 + 9.0*p11 + 9.0*p22 - p33)/16.0 derived from Nevilles's algorithm
    //interpolate real part
    return (-((double*)&(this->Xt.at((lower_hist+0)%this->hist_length)))[index] + 9.0*((double*)&(this->Xt.at((lower_hist+1)%this->hist_length)))[index] + 9.0*((double*)&(this->Xt.at((lower_hist+2)%this->hist_length)))[index] - ((double*)&(this->Xt.at((lower_hist+3)%this->hist_length)))[index])/16.0;
  }
  
  std::complex<double> at_IntPol_3rd_C_hlfstp(double tau, int index){
  //returns a complex<double> of the indexed variable of vars which is interpolated to third order using Neville's algorithm at the time t0-tau
    int lower_hist = (this->t0 + 2*this->hist_length - (int)floor((tau-this->tau_offset)/this->dt) - 2) % this->hist_length;
    //interpolate according to p03 = (-p00 + 9.0*p11 + 9.0*p22 - p33)/16.0 derived from Nevilles's algorithm
    //interpolate real part
    const double p03R = (-((double*)&(this->Xt[(lower_hist+0)%this->hist_length]))[index] + 9.0*((double*)&(this->Xt[(lower_hist+1)%this->hist_length]))[index] + 9.0*((double*)&(this->Xt[(lower_hist+2)%this->hist_length]))[index] - ((double*)&(this->Xt[(lower_hist+3)%this->hist_length]))[index])/16.0;
    //interpolate imaginary part
    const double p03C = (-((double*)&(this->Xt[(lower_hist+0)%this->hist_length]))[index+1] + 9.0*((double*)&(this->Xt[(lower_hist+1)%this->hist_length]))[index+1] + 9.0*((double*)&(this->Xt[(lower_hist+2)%this->hist_length]))[index+1] - ((double*)&(this->Xt[(lower_hist+3)%this->hist_length]))[index+1])/16.0;
    
    return std::complex<double>(p03R,p03C);
  }

      
};


class vars_vec_wdX : public vars_vec{

private:
  
public:
  std::vector<vars> dXtdt;
  
  //constructor stuff
  vars_vec_wdX(double tau_max, double dt) : vars_vec(tau_max,dt) {
    dXtdt.resize(hist_length);
  }
  
  using vars_vec::vars_vec;
  
  void setToCnst(double s){
    //sets all variables of all vars to the input s
    for (auto & v : this->Xt) {
      v.setTo(s);
    }
    for (auto & v : this->dXtdt) {
      v.setTo(s);
    }
  }
  
  vars * dX_at_t0_rPtr(){
    //returns pointer to vars at current time t0
    return &(this->dXtdt[this->t0]);
  }
  
  vars * dX_at_nxt_rPtr(){
    //returns pointer to the newest (next) vars, which previously was the oldest vars in the history
    return &(this->dXtdt[(this->t0 + 1)%this->hist_length]);
  }

  vars * dX_at_Rnd_rPtr(double tau){
    //return pointer to vars at time t0-tau where tau/dt is rounded to the closest integer
    if(tau > this->tau_max){
      std::cout << "Error: requested delay tau is greater than tau_max" << std::endl;
      return nullptr;
    }
    else{
      return &(this->dXtdt[(this->t0 + this->hist_length - (int)round((tau-this->tau_offset)/this->dt)) % this->hist_length]);
    }
  }
  
  void save(std::string filename){
    std::ofstream binout_Xhist (filename, std::ios::out | std::ios::binary);
    
    for(std::size_t i=0; i<sizeof(int); i++) binout_Xhist << ((unsigned char*)(&(hist_length)))[i];
    for(std::size_t i=0; i<sizeof(int); i++) binout_Xhist << ((unsigned char*)(&(t0)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(dt)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(tau_max)))[i];
    for(std::size_t i=0; i<sizeof(double); i++) binout_Xhist << ((unsigned char*)(&(tn)))[i];
    
    for(unsigned int k = 0; k < hist_length; k++){
      for(unsigned int l = 0; l < sizeof(vars); l++){
        binout_Xhist << ((unsigned char*)(&this->Xt[k]))[l];
      }
    }
    
    for(unsigned int k = 0; k < hist_length; k++){
      for(unsigned int l = 0; l < sizeof(vars); l++){
        binout_Xhist << ((unsigned char*)(&this->dXtdt[k]))[l];
      }
    }
    
    binout_Xhist.close();
  }
  
  void load(std::string filename){
    std::ifstream binin_Xhist (filename, std::ios::in | std::ios::binary);
    if(binin_Xhist.good()){
    
      for(std::size_t i=0; i<sizeof(int); i++) ((unsigned char*)(&(hist_length)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(int); i++) ((unsigned char*)(&(t0)))[i] = binin_Xhist.get();
      this->Xt.resize(this->hist_length);
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(dt)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(tau_max)))[i] = binin_Xhist.get();
      for(std::size_t i=0; i<sizeof(double); i++) ((unsigned char*)(&(tn)))[i] = binin_Xhist.get();
      
      for(unsigned int k = 0; k < hist_length; k++){
        for(unsigned int l = 0; l < sizeof(vars); l++){
          ((unsigned char*)(&this->Xt[k]))[l] = binin_Xhist.get();
        }
      }
      
      for(unsigned int k = 0; k < hist_length; k++){
        for(unsigned int l = 0; l < sizeof(vars); l++){
          ((unsigned char*)(&this->dXtdt[k]))[l] = binin_Xhist.get();
        }
      }

    }
    else std::cout << "Error: History file not found" << std::endl;
    
    binin_Xhist.close();
  }
  
  double at_R_cubicHermite(double tau, int index){
    int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
    int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
    double rel_dt = (tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt);
    
    const double a0 = 2.0 * rel_dt * rel_dt * rel_dt - 3.0 * rel_dt * rel_dt + 1.0;
    const double a1 = -2.0 * rel_dt * rel_dt * rel_dt + 3.0 * rel_dt * rel_dt;
    const double b0 = rel_dt * rel_dt * rel_dt - 2.0 * rel_dt * rel_dt + rel_dt;
    const double b1 = rel_dt * rel_dt * rel_dt - rel_dt * rel_dt;
    
    const double f0 = ((double*)&(this->Xt[lower_hist]))[index];
    const double f1 = ((double*)&(this->Xt[upper_hist]))[index];
    const double df0 = ((double*)&(this->dXtdt[lower_hist]))[index];
    const double df1 = ((double*)&(this->dXtdt[upper_hist]))[index];
    
    return f0*a0 + f1*a1 - this->dt*(df0*b0 + df1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.
  }
  

  std::complex<double> at_C_cubicHermite(double tau, int index){
    int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
    int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
    double rel_dt = (tau-this->tau_offset)/this->dt - floor((tau-this->tau_offset)/this->dt);
    
    const double a0 = 2.0 * rel_dt * rel_dt * rel_dt - 3.0 * rel_dt * rel_dt + 1.0;
    const double a1 = -2.0 * rel_dt * rel_dt * rel_dt + 3.0 * rel_dt * rel_dt;
    const double b0 = rel_dt * rel_dt * rel_dt - 2.0 * rel_dt * rel_dt + rel_dt;
    const double b1 = rel_dt * rel_dt * rel_dt - rel_dt * rel_dt;
    
    const double f0 = ((double*)&(this->Xt[lower_hist]))[index];
    const double f1 = ((double*)&(this->Xt[upper_hist]))[index];
    const double df0 = ((double*)&(this->dXtdt[lower_hist]))[index];
    const double df1 = ((double*)&(this->dXtdt[upper_hist]))[index];
    
    const double real = f0*a0 + f1*a1 - this->dt*(df0*b0 + df1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.
    
    //imag part
    const double if0 = ((double*)&(this->Xt[lower_hist]))[index+1];
    const double if1 = ((double*)&(this->Xt[upper_hist]))[index+1];
    const double idf0 = ((double*)&(this->dXtdt[lower_hist]))[index+1];
    const double idf1 = ((double*)&(this->dXtdt[upper_hist]))[index+1];
    
    const double imag = if0*a0 + if1*a1 - this->dt*(idf0*b0 + idf1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.
    
    return std::complex<double>(real,imag);
  }

  
  //////////////////////////////////////
  // fast cubic hermite interpol with halfsteps
  //////////////////////////////////////
  
  double at_R_fCH(double tau, int index){
    if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
      return at_IntPol_CH_R_hlfstp(tau, index);
    }
    else{
      return ((double*)this->at_Rnd_rPtr(tau))[index];
    }
  }
  
  std::complex<double> at_C_fCH(double tau, int index){
    if(this->tau_offset > 0.3*dt && this->tau_offset < 0.7*dt){
      return at_IntPol_CH_C_hlfstp(tau, index);
    }
    else{
      int rind = (this->t0 + this->hist_length - (int)round((tau-this->tau_offset)/this->dt)) % this->hist_length;
      return std::complex<double>(((double*)(&(this->Xt[rind])))[index],((double*)(&(this->Xt[rind])))[index+1]);
    }
  }
  
  //////////////////////////////////////
  // only to be used to calc cubic hermite halfsteps
  //////////////////////////////////////
  
  double at_IntPol_CH_R_hlfstp(double tau, int index){
    //cubic hermite interpolation at rel_dt = 0.5
    
    const int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
    const int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
    
    const double a0 = 0.5;
    const double a1 = 0.5;
    const double b0 = 0.125;
    const double b1 = -0.125;
    
    const double f0 = ((double*)&(this->Xt[lower_hist]))[index];
    const double f1 = ((double*)&(this->Xt[upper_hist]))[index];
    const double df0 = ((double*)&(this->dXtdt[lower_hist]))[index];
    const double df1 = ((double*)&(this->dXtdt[upper_hist]))[index];
    
    return f0*a0 + f1*a1 - this->dt*(df0*b0 + df1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.

  }
  
  std::complex<double> at_IntPol_CH_C_hlfstp(double tau, int index){
    //cubic hermite interpolation at rel_dt = 0.5
    
    const int lower_hist = (this->t0 + this->hist_length - (int)floor((tau-this->tau_offset)/this->dt)) % this->hist_length;
    const int upper_hist = (lower_hist + this->hist_length - 1)%this->hist_length;
    
    const double a0 = 0.5;
    const double a1 = 0.5;
    const double b0 = 0.125;
    const double b1 = -0.125;
    
    const double f0 = ((double*)&(this->Xt[lower_hist]))[index];
    const double f1 = ((double*)&(this->Xt[upper_hist]))[index];
    const double df0 = ((double*)&(this->dXtdt[lower_hist]))[index];
    const double df1 = ((double*)&(this->dXtdt[upper_hist]))[index];
    
    const double real = f0*a0 + f1*a1 - this->dt*(df0*b0 + df1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.
    
    //imag part
    const double if0 = ((double*)&(this->Xt[lower_hist]))[index+1];
    const double if1 = ((double*)&(this->Xt[upper_hist]))[index+1];
    const double idf0 = ((double*)&(this->dXtdt[lower_hist]))[index+1];
    const double idf1 = ((double*)&(this->dXtdt[upper_hist]))[index+1];
    
    const double imag = if0*a0 + if1*a1 - this->dt*(idf0*b0 + idf1*b1); //!!!!!! minus at the df terms as time in the history is backwards -> derivatives change sign.
    
    return std::complex<double>(real,imag);

  }
  
};







