// #define stoch_spon
#define det_spon
// #define kick

#include "incl/global_v2.cpp"
#include "incl/get_params.cc"
#include "incl/get_specs.cc"
#include "incl/get_correlations.cc"
#include "TSQDL_parameters.cpp"
#include "TSQDL_vars.cpp"
#include "incl/vars_vec_v8.cpp"
#include "incl/DDEintegrator_v4.cpp"
#include "incl/evaluateTS_v2.cpp"
#include "TSQDL_eqs.cpp"


using namespace std;


int main(int argc, char* argv[]){
    
  //paramters for the QD model in LIN14units
  parameters p;
  p.gGS = katana::getCmdOption(argv, argv+argc, "-gGS" , 230.0);			//ground state gain
  p.gES = katana::getCmdOption(argv, argv+argc, "-gES" , 0.0); 			//excited state gain
  p.kappaGS = katana::getCmdOption(argv, argv+argc, "-kappaGS" , 50.0);			//optical losses for the GS mode
  p.kappaES = katana::getCmdOption(argv, argv+argc, "-kappaES" , 50.0);			//optical losses for the ES mode
  p.beta = katana::getCmdOption(argv, argv+argc, "-beta" , 2.2E-3);			//chance to spontaneously emit into lasing mode
  p.NQD = katana::getCmdOption(argv, argv+argc, "-NQD" , 1.0);				//area density of QDs in
  p.ZQD = katana::getCmdOption(argv, argv+argc, "-ZQD" , 3.0E7);			//number of QDs
  p.fact = katana::getCmdOption(argv, argv+argc, "-fact" , 0.5);			//fraction of active QDs
  p.WGS = katana::getCmdOption(argv, argv+argc, "-WGS" , 0.44);				//Einstein coefficient for spontaneous emission for the GS
  p.WES = katana::getCmdOption(argv, argv+argc, "-WES" , 0.55);				//Einstein coefficient for spontaneous emission for the ES
  p.RWloss = katana::getCmdOption(argv, argv+argc, "-RWloss" , 0.54);			//loss rate for the quantum wells
  p.etaGS = katana::getCmdOption(argv, argv+argc, "-etaGS" , 10.055E-6*sqrt(10));	//electric field conversion factor for the GS
  p.etaES = katana::getCmdOption(argv, argv+argc, "-etaES" , 10.416E-6*sqrt(10));	//electric field conversion factor for the ES
  p.GS_del_om_ES = katana::getCmdOption(argv, argv+argc, "-GS_del_om_ES" , 125.0);	//fits from LIN14	
  p.GS_del_om_QW_e = katana::getCmdOption(argv, argv+argc, "-GS_del_om_QW_e" , 11.0);	//fits from LIN14
  p.GS_del_om_QW_h = katana::getCmdOption(argv, argv+argc, "-GS_del_om_QW_h" , 5.5);	//fits from LIN14
  p.ES_del_om_GS = katana::getCmdOption(argv, argv+argc, "-ES_del_om_GS" , -62.5);	//half of what the GS sees due to half the degeneracy		
  p.ES_del_om_QW_e = katana::getCmdOption(argv, argv+argc, "-ES_del_om_QW_e" , 15.0);	//slightly more than the GS
  p.ES_del_om_QW_h = katana::getCmdOption(argv, argv+argc, "-ES_del_om_QW_h" , 7.0);	//slightly more than the GS
  p.gamma_GS = katana::getCmdOption(argv, argv+argc, "-gamma_GS" , 1.0);		//scales GS amplitude-phase coupling
  p.gamma_ES = katana::getCmdOption(argv, argv+argc, "-gamma_ES" , 0.0);		//scales ES amplitude-phase coupling
  
  p.SScale = katana::getCmdOption(argv, argv+argc, "-SScale" , 1.0);
  
  //controll parameters
  p.J = katana::getCmdOption(argv, argv+argc, "-J" , 28.0);
  p.tau_GS = katana::getCmdOption(argv, argv+argc, "-tau_GS" , 3.6);
  p.tau_ES = katana::getCmdOption(argv, argv+argc, "-tau_ES" , 0.15);
  p.K_GS = katana::getCmdOption(argv, argv+argc, "-K_GS" , 0.0);
  p.K_ES = katana::getCmdOption(argv, argv+argc, "-K_ES" , 0.0);
  
//   p.C_GS = katana::getCmdOption(argv, argv+argc, "-C_GS" , 4.602);
//   p.img_offset = 184.423; // old
  
//   p.C_GS = katana::getCmdOption(argv, argv+argc, "-C_GS" , 0.8495234641020676153735);
//   p.img_offset = 0; // old
  
  p.C_GS = katana::getCmdOption(argv, argv+argc, "-C_GS" , 3.14678119); //trial and error
//   p.C_GS = katana::getCmdOption(argv, argv+argc, "-C_GS" , 3.14679257); //calculated
  p.img_offset = 236.012284903; //J=28
  
  p.C_ES = katana::getCmdOption(argv, argv+argc, "-C_ES" , 0);

  //prevent GS or ES from going to AND staying at 0
  p.fakeSpE = katana::getCmdOption(argv, argv+argc, "-fakeSpe" , 1.0E-50);
  p.kick_time = katana::getCmdOption(argv, argv+argc, "-kick_time" , 10.0);
  p.kick_strength = katana::getCmdOption(argv, argv+argc, "-kick_strength" , 1.0E-6);
  
	
  //solver
  DDEintegrator DDEsolver;
  //solver parameters
  unsigned long long int tn = 0;
  double intTime = katana::getCmdOption(argv, argv+argc, "-intTime" , 10.0);
  double dt = katana::getCmdOption(argv, argv+argc, "-dt" , 0.001);
  double outTime = katana::getCmdOption(argv, argv+argc, "-outTime" , 10.0);
  if(outTime < 0.0) outTime = intTime;
  unsigned long long int outTime_ntn = (unsigned long long int)(outTime/dt);
  
  
  //dynamical variables with delay
//   vars_vec Xhist(std::max(p.tau_ES, p.tau_GS),dt);
  vars_vec_wdX Xhist(std::max(p.tau_ES, p.tau_GS),dt);
  
  
  //////////////////////////////////////////
  //noise
  //////////////////////////////////////////
  
  katana::Seed s;
  s.set_seed();
  
  //noise
  p.SqrtNoiseStr_GS = sqrt( katana::getCmdOption(argv, argv+argc, "-noiseStr_GS" , 1.0) * (p.beta * p.ZQD * p.fact * p.WGS * p.etaGS * p.etaGS) / dt );
  p.SqrtNoiseStr_ES = sqrt( katana::getCmdOption(argv, argv+argc, "-noiseStr_ES" , 1.0) * (p.beta * p.ZQD * p.fact * p.WES * p.etaES * p.etaES) / dt );
  
//   std::function<void (parameters*)> noise = noise_empty; //default = no noise
//     
//   //gw noise
//   if(katana::getCmdOption_bool(argv, argv+argc, "-noise" , false)){
//     noise = [&](parameters* p){ 
//       p->noise_GS = p->SqrtNoiseStr_GS*katana::gwNoise();
//       p->noise_ES = p->SqrtNoiseStr_ES*katana::gwNoise();
//     };
//   }
  
  std::function<void (parameters*)> noise = [&](parameters* p){ 
      p->noise_GS = p->SqrtNoiseStr_GS*katana::gwNoise();
      p->noise_ES = p->SqrtNoiseStr_ES*katana::gwNoise();
  };
    
  //gw noise
  if(katana::getCmdOption_bool(argv, argv+argc, "-nonoise" , false)) noise = noise_empty;
  
  
  //////////////////////////////////////////
  //outputfunctions
  //////////////////////////////////////////
  
  //for timeseries output
  outputfile.precision(10); //outputfile must be declared in global.hpp
  
  //empty output function
  auto empty = [&](vars* X, double t, double TMax){};
  
  //output GS and ES
  auto outputLRToFile = [&](vars* X, double t, double TMax){
    //column labels
    if(t==0){
      outputfile << "#t(1) GS(2) ES(3)" << std::endl;
    }
    outputfile << t << "\t" << norm(X->EGS) << "\t" <<  norm(X->EES) << std::endl;
  };
  
  //output timeseries of time and all dynamical variables -> very slow and might take lots  of memory
  auto outputAllToFile = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
    int k = 0;
    int kmax = 1;
    if(tn >= tn_final - outTime_ntn && k == 0){
      outputfile << tn*dt << "\t";
      outputfile << norm(X->EGS) << "\t" << norm(X->EES) << "\t" << X->rhoGSe << "\t" << X->rhoGSh << "\t" << X->rhoESe << "\t" << X->rhoESh << "\t" << X->rhoGSei << "\t" << X->rhoGShi << "\t" << X->rhoESei << "\t" << X->rhoEShi << "\t" << X->we << "\t" << X->wh << "\t" << real(X->EGS) << "\t" << imag(X->EGS);
      outputfile << std::endl;
    }
    k = (k+1)%kmax;
  };
  
  //output to vectors for TS analysis
  std::vector<double> outVec_GS, outVec_ES;
  outVec_GS.reserve((int)(outTime/dt)), outVec_ES.reserve((int)(outTime/dt));
  
  std::vector<std::complex<double>> outVecComplex_GS, outVecComplex_ES;
  outVecComplex_GS.reserve((int)(outTime/dt)), outVecComplex_ES.reserve((int)(outTime/dt));
  
  std::vector<double> GSinv; GSinv.reserve((int)(outTime/dt));
  std::vector<double> ESinv; ESinv.reserve((int)(outTime/dt));
  
  auto outputToVector = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn){
      outVec_GS.push_back(norm(X->EGS)*p.kappaGS*0.7242); //output in Watts
      outVec_ES.push_back(norm(X->EES)*p.kappaES*0.7242);
      outVecComplex_GS.push_back(X->EGS);
//       outVecComplex_ES.push_back(X->EES);
      
      GSinv.push_back(X->rhoGSe + X->rhoGSh - 1.0);
      ESinv.push_back(X->rhoESe + X->rhoESh - 1.0);
    }
  };
  
  
  //////////////////////////////////////////
  //for sweeps
  //////////////////////////////////////////
  
  std::ofstream out_sweep_State, out_sweep_Ext_GS_max, out_sweep_Ext_GS_min, out_sweep_Ext_ES_max, out_sweep_Ext_ES_min, out_sweep_TS;
  std::ofstream out_sweep_AllExt_GS_max, out_sweep_AllExt_GS_min, out_sweep_AllExt_ES_max, out_sweep_AllExt_ES_min;
  out_sweep_State.precision(12), out_sweep_TS.precision(12), out_sweep_Ext_GS_max.precision(12), out_sweep_Ext_GS_min.precision(12), out_sweep_Ext_ES_max.precision(12), out_sweep_Ext_ES_min.precision(12);
  out_sweep_AllExt_GS_max.precision(12), out_sweep_AllExt_GS_min.precision(12), out_sweep_AllExt_ES_max.precision(12), out_sweep_AllExt_ES_min.precision(12);
  std::ofstream out_sweep_Ext_freq_max, out_sweep_Ext_freq_min;
  out_sweep_Ext_freq_max.precision(12), out_sweep_Ext_freq_min.precision(12);
  evaluateTS eval_sweep_GS;
  evaluateTS eval_sweep_ES;
  evaluateTS eval_sweep_freq;
  
  double sweep_doubleCountTol = katana::getCmdOption(argv, argv+argc, "-SDoubleCountTol" , 1E-2);
  eval_sweep_GS.doubleCountTol = sweep_doubleCountTol;
  eval_sweep_ES.doubleCountTol = sweep_doubleCountTol;
  eval_sweep_freq.doubleCountTol = sweep_doubleCountTol;
  
  double* PtrToSweepPar;
  std::string str_sweep_parameter;
  double* PtrToSecPar;
  std::string str_sweep_sec_parameter;
  
  double sweep_IntTime = katana::getCmdOption(argv, argv+argc, "-sIntTime" , 1000.0);
  double sweep_OutTime = katana::getCmdOption(argv, argv+argc, "-sOutTime" , 100.0);
  int nsweep_steps = katana::getCmdOption(argv, argv+argc, "-sSteps" , 100);
  double sweep_start = katana::getCmdOption(argv, argv+argc, "-sStart" , 0);
  double sweep_end = katana::getCmdOption(argv, argv+argc, "-sEnd" , 1);
  double sweep_step = (sweep_end-sweep_start)/(double)nsweep_steps;
  
  std::vector<double> sweep_Phase, sweep_Frequency, sweep_rollingMeanFreq;
  sweep_Phase.reserve((int)(outTime/dt)), sweep_Frequency.reserve((int)(outTime/dt)), sweep_rollingMeanFreq.reserve((int)(outTime/dt));
  
  std::string str_sweep_data = "data/";
  std::string str_sweep_updown = "up";
  
  std::string str_state;
  std::string str_extrema_GS;
  std::string str_extrema_ES;
  std::string str_extrema_freq;
  std::string str_powerSpec;
  std::string str_opticalSpec;
  std::string str_AC;
  std::string str_TS;
  std::string str_bin;
  
  bool bool_wdownsweep = katana::getCmdOption_bool(argv, argv+argc, "-wdownSweep" , false);
  bool bool_powerSpec = katana::getCmdOption_bool(argv, argv+argc, "-wpowerSpec" , false);
  bool bool_opticalSpec = katana::getCmdOption_bool(argv, argv+argc, "-wopticalSpec" , false);
  bool bool_AC = katana::getCmdOption_bool(argv, argv+argc, "-wAC" , false);
  bool bool_TS = katana::getCmdOption_bool(argv, argv+argc, "-wTS" , false);
  bool bool_Extrema = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema" , false);
  bool bool_AllExtrema = katana::getCmdOption_bool(argv, argv+argc, "-wAllExtrema" , false);
  bool bool_bin = katana::getCmdOption_bool(argv, argv+argc, "-wBin" , false);
  

  
  auto sweep = [&](){
    
    str_state = str_sweep_data + str_sweep_updown + "_state_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_extrema_GS = str_sweep_data + str_sweep_updown + "_extrema_GS_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_extrema_ES = str_sweep_data + str_sweep_updown + "_extrema_ES_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_extrema_freq = str_sweep_data + str_sweep_updown + "_extrema_freq_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_powerSpec = str_sweep_data + "/powerSpec/powerSpec_"+ str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_opticalSpec = str_sweep_data + "/opticalSpec/opticalSpec_"+ str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_AC = str_sweep_data + "AC/AC_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_TS = str_sweep_data + "TS/TS_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_bin = str_sweep_data + "bin/bin_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    
    if(bool_Extrema){
      out_sweep_Ext_GS_max.open(str_extrema_GS + "_max");
      out_sweep_Ext_GS_min.open(str_extrema_GS + "_min");
      out_sweep_Ext_ES_max.open(str_extrema_ES + "_max");
      out_sweep_Ext_ES_min.open(str_extrema_ES + "_min");
      out_sweep_Ext_freq_max.open(str_extrema_freq + "_max");
      out_sweep_Ext_freq_min.open(str_extrema_freq + "_min");
    }
    if(bool_AllExtrema){
      out_sweep_AllExt_GS_max.open(str_extrema_GS + "_All_max");
      out_sweep_AllExt_GS_min.open(str_extrema_GS + "_All_min");
      out_sweep_AllExt_ES_max.open(str_extrema_ES + "_All_max");
      out_sweep_AllExt_ES_min.open(str_extrema_ES + "_All_min");
    }
    out_sweep_State.open(str_state);
    
    outTime = sweep_OutTime;
    *PtrToSweepPar = sweep_start;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-9);
    }
    
    for(int k = 0; k <= nsweep_steps; k++){
      sm::clck = clock();
      tn=0;
//       Xhist.setToCnst(1E-6);
      outVec_GS.resize(0), outVec_ES.resize(0), outVecComplex_GS.resize(0), outVecComplex_ES.resize(0), GSinv.resize(0), ESinv.resize(0);
      DDEsolver.DDE_RK4(TSQDL, algs_empty, noise, &Xhist, &p, &tn, sweep_IntTime, dt, outputToVector);
      eval_sweep_GS.newTS(&outVec_GS, dt);
      eval_sweep_ES.newTS(&outVec_ES, dt);
      
      //analize optical frequency
      
      sweep_Phase.resize(0), sweep_Frequency.resize(0), sweep_rollingMeanFreq.resize(0);
      double smfreq =  eval_sweep_GS.opticalFreq(outVecComplex_GS, sweep_Phase, sweep_Frequency, dt);
      int sk_tau = p.tau_GS / dt;
      
      for(int k = sk_tau; k < sweep_Phase.size()-2; k++){
        sweep_rollingMeanFreq.push_back( (sweep_Phase[k] - sweep_Phase[k-sk_tau])/p.tau_GS );
      }
      eval_sweep_freq.newTS(&sweep_rollingMeanFreq,dt);
      
//       out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//       for(int k = 0; k < sweep_rollingMeanFreq.size(); k++) out_sweep_TS << dt*k << "\t" << sweep_rollingMeanFreq[k] << std::endl;
//       out_sweep_TS.close();
 
      cout << str_sweep_parameter << ": " << *PtrToSweepPar;
//       out_sweep_State << p.J << "\t" << p.K_GS << "\t" << p.tau_GS << "\t" << p.C_GS << "\t";
      out_sweep_State << *PtrToSweepPar << "\t" << *PtrToSecPar << "\t"; //2
      out_sweep_State << eval_sweep_GS.average << "\t" << eval_sweep_GS.greatestMax << "\t" << eval_sweep_GS.smallestMin << "\t"; //5
      out_sweep_State << eval_sweep_GS.numberOfMax() << "\t" << eval_sweep_GS.numberOfMin() << "\t"; //7
      
      out_sweep_State << eval_sweep_ES.average << "\t" << eval_sweep_ES.greatestMax << "\t" << eval_sweep_ES.smallestMin << "\t"; //10
      out_sweep_State << eval_sweep_ES.numberOfMax() << "\t" << eval_sweep_ES.numberOfMin() << "\t"; //12
      out_sweep_State << smfreq << "\t"; //13
//       out_sweep_State << eval_sweep_ES.opticalFreq(outVecComplex_ES, sweep_Phase, sweep_Frequency, dt) << "\t";
      out_sweep_State << eval_sweep_freq.average << "\t" << eval_sweep_freq.numberOfMax() << "\t" << eval_sweep_freq.numberOfMin() << "\t" << eval_sweep_freq.loopsRatio() << "\t"; //17
      out_sweep_State << endl;
      
      if(bool_Extrema){
        //GS
        out_sweep_Ext_GS_max << *PtrToSweepPar << "\t" << eval_sweep_GS.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_GS.uniqueMax.size(); k++) out_sweep_Ext_GS_max << *PtrToSweepPar << "\t" << eval_sweep_GS.uniqueMax[k] << std::endl;      
        out_sweep_Ext_GS_min << *PtrToSweepPar << "\t" <<  eval_sweep_GS.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep_GS.uniqueMin.size(); k++) out_sweep_Ext_GS_min << *PtrToSweepPar << "\t" << eval_sweep_GS.uniqueMin[k] << std::endl;
        //ES
        out_sweep_Ext_ES_max << *PtrToSweepPar << "\t" << eval_sweep_ES.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_ES.uniqueMax.size(); k++) out_sweep_Ext_ES_max << *PtrToSweepPar << "\t" << eval_sweep_ES.uniqueMax[k] << std::endl;      
        out_sweep_Ext_ES_min << *PtrToSweepPar << "\t" <<  eval_sweep_ES.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep_ES.uniqueMin.size(); k++) out_sweep_Ext_ES_min << *PtrToSweepPar << "\t" << eval_sweep_ES.uniqueMin[k] << std::endl;
        //Frequency
        out_sweep_Ext_freq_max<< *PtrToSweepPar << "\t" << eval_sweep_freq.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_freq.uniqueMax.size(); k++) out_sweep_Ext_freq_max << *PtrToSweepPar << "\t" << eval_sweep_freq.uniqueMax[k] << std::endl;      
        for(int k = 0; k < eval_sweep_freq.uniqueMin.size(); k++) out_sweep_Ext_freq_min << *PtrToSweepPar << "\t" << eval_sweep_freq.uniqueMin[k] << std::endl;
        out_sweep_Ext_freq_min << *PtrToSweepPar << "\t" <<  eval_sweep_freq.smallestMin << std::endl;
      }
      if(bool_AllExtrema){
        //GS
        out_sweep_AllExt_GS_max << *PtrToSweepPar << "\t" << eval_sweep_GS.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_GS.maxima.size(); k++) out_sweep_AllExt_GS_max << *PtrToSweepPar << "\t" << eval_sweep_GS.maxima[k] << std::endl;      
        out_sweep_AllExt_GS_min << *PtrToSweepPar << "\t" <<  eval_sweep_GS.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep_GS.minima.size(); k++) out_sweep_AllExt_GS_min << *PtrToSweepPar << "\t" << eval_sweep_GS.minima[k] << std::endl;
        //ES
        out_sweep_AllExt_ES_max << *PtrToSweepPar << "\t" << eval_sweep_ES.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_ES.maxima.size(); k++) out_sweep_AllExt_ES_max << *PtrToSweepPar << "\t" << eval_sweep_ES.maxima[k] << std::endl;      
        out_sweep_AllExt_ES_min << *PtrToSweepPar << "\t" <<  eval_sweep_ES.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep_ES.minima.size(); k++) out_sweep_AllExt_ES_min << *PtrToSweepPar << "\t" << eval_sweep_ES.minima[k] << std::endl;
      }
      if(bool_powerSpec){
        std::vector<double> powerSpec;
        std::vector<double> WindowedData;
        //GS
        sm::window_Hann(outVec_GS, WindowedData);
        katana::get_power_spec(WindowedData, powerSpec);
        std::ostringstream sw_num;
        sw_num << std::fixed << std::setprecision(12);
        sw_num << *PtrToSweepPar;
        katana::sm_dump_power_spec(powerSpec, dt, 10.0, str_powerSpec+"_GS_Hann_"+str_sweep_parameter+"_"+sw_num.str());
//         katana::sm_dump_power_spec(powerSpec, dt, 1.0, str_powerSpec+"_GS_windowed_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//         //no window
//         katana::get_power_spec(outVec_GS, powerSpec);
//         katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, str_powerSpec+"_GS_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//         //ES
//         sm::window_Hann(outVec_ES, WindowedData);
//         katana::get_power_spec(WindowedData, powerSpec);
//         katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, str_powerSpec+"_ES_windowed_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//         //no window
//         katana::get_power_spec(outVec_ES, powerSpec);
//         katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, str_powerSpec+"_ES_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
      }
      if(bool_opticalSpec){
        std::vector<double> opticalSpec;
        std::ostringstream sw_num;
        sw_num << std::fixed << std::setprecision(12);
        sw_num << *PtrToSweepPar;
        
        std::vector<std::complex<double>> WindowedData;
        sm::window_Hann(outVecComplex_GS, WindowedData);
        katana::get_optical_spec(WindowedData, opticalSpec);
        katana::sm_dump_optical_spec(opticalSpec, dt, 4, str_opticalSpec+"_GS_Hann_"+str_sweep_parameter+"_"+sw_num.str());
//         //no window
//         katana::get_optical_spec(outVecComplex_GS, opticalSpec);
//         katana::dump_optical_spec(opticalSpec, dt, opticalSpec.size()/10, str_opticalSpec+"_GS_"+str_sweep_parameter+"_"+sw_num.str());
      }
      if(bool_TS){
        out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
        for(int k = 0; k < outVec_GS.size()/10; k++) out_sweep_TS << dt*k << "\t" << outVec_GS[k] << "\t" << outVec_ES[k] << std::endl;
        out_sweep_TS.close();
      }
      if(bool_AC){
        std::vector<double> AC;
        //GS
        katana::get_autocorr_norm_out(outVecComplex_GS, AC);
        katana::print_correlation(AC, dt, str_AC+"_GS_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar), (int)(1000.0/dt), 1, "ps");
        //ES
        katana::get_autocorr_norm_out(outVecComplex_ES, AC);
        katana::print_correlation(AC, dt, str_AC+"_ES_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar), (int)(1000.0/dt), 1, "ps");
      }
      if(bool_bin){
        Xhist.tn = tn;
        Xhist.save(str_bin+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
        cout << str_bin+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar) << endl;
      }
      
      cout << ", integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
      *PtrToSweepPar += sweep_step;
    }
    if(bool_Extrema){
      out_sweep_Ext_GS_max.close();
      out_sweep_Ext_GS_min.close();
      out_sweep_Ext_ES_max.close();
      out_sweep_Ext_ES_min.close();
      out_sweep_Ext_freq_max.close();
      out_sweep_Ext_freq_min.close();
    }
    if(bool_AllExtrema){
      out_sweep_AllExt_GS_max.close();
      out_sweep_AllExt_GS_min.close();
      out_sweep_AllExt_ES_max.close();
      out_sweep_AllExt_ES_min.close();
    }
    out_sweep_State.close();
  };
 
  
  
  //timeseries to file
  if(katana::getCmdOption_bool(argv, argv+argc, "-allVarsTS" , false)){
    outputfile.open("data/allVarsTS");
    tn=0;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-9);
    }
    
    sm::clck = clock();
    DDEsolver.DDE_RK4(TSQDL, algs_empty, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
//     DDEsolver.DDE_euler(MLL_derivs, algs_empty, noise, &Xhist, &p, &t, T, dt, outputAllToFile);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    outputfile.close();
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
  }
  
  
  
  //compute and analize time series
  if(katana::getCmdOption_bool(argv, argv+argc, "-TS" , false)){
    
    //evaluateTS objects
    evaluateTS eval_GS;
    evaluateTS eval_ES;
    evaluateTS eval_avfreq;
    outVec_GS.resize(0), outVec_ES.resize(0), outVecComplex_GS.resize(0), outVecComplex_ES.resize(0);
    
    std::vector<double> Phase, Frequency;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(0.001);
    }

    sm::clck = clock();
//     DDEsolver.DDE_euler(TSQDL, algs_empty, noise, &Xhist, &p, &t, T, dt, outputToVector);
    DDEsolver.DDE_RK4(TSQDL, algs_empty, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    //load new TS
    eval_GS.newTS(&outVec_GS, dt);
    eval_ES.newTS(&outVec_ES, dt);
    
    double avFreq = eval_GS.opticalFreq(outVecComplex_GS, Phase, Frequency, dt);
    int k_tau = p.tau_GS / dt;
    std::vector<double> avfreq;
    avfreq.reserve(Phase.size());
    for(int k = k_tau; k < Phase.size()-2; k++){
      avfreq.push_back( (Phase[k] - Phase[k-k_tau])/p.tau_GS );
    }
    eval_avfreq.newTS(&avfreq,dt);
    
    cout << "----results----" << endl;
    cout <<"GS average: " << eval_GS.average << endl;
    cout <<"ES average: " << eval_ES.average << endl;
    cout <<"GS greates Max: " << eval_GS.greatestMax << ", GS smallest Min: " << eval_GS.smallestMin << endl;
    cout <<"ES greates Max: " << eval_ES.greatestMax << ", GS smallest Min: " << eval_ES.smallestMin << endl;
    cout <<"GS number of Max: " << eval_GS.numberOfMax() <<", GS number of Min: " << eval_GS.numberOfMin() << endl;
    cout <<"ES number of Max: " << eval_ES.numberOfMax() <<", ES number of Min: " << eval_ES.numberOfMin() << endl;
    cout <<"GS period: " << eval_GS.period << endl;
    cout <<"ES period: " << eval_ES.period << endl;
    cout <<"GS average optical frequency: " << avFreq << endl;
    cout <<"Rolling mean frequency: mean: " << eval_avfreq.average << " greatest max: " << eval_avfreq.greatestMax << " smallest min: " << eval_avfreq.smallestMin << " number of max: " << eval_avfreq.numberOfMax() << " number of min: " << eval_avfreq.numberOfMin() << endl;
    cout <<"loop ratio: " <<  eval_avfreq.loopsRatio() << endl;
    cout << "----^^^^^^----" << endl;
//     cout << std::setprecision(12) << p.gamma_GS * ( p.GS_del_om_ES * (p.fact *(Xhist.at_t0_rPtr()->rhoESe + Xhist.at_t0_rPtr()->rhoESh) + (1-p.fact)*(Xhist.at_t0_rPtr()->rhoESei + Xhist.at_t0_rPtr()->rhoEShi)) + p.GS_del_om_QW_e * Xhist.at_t0_rPtr()->we + p.GS_del_om_QW_h * Xhist.at_t0_rPtr()->wh) << endl;
    
    outputfile.open("data/out_TS");
    for(int k = 0; k < outVec_GS.size()-2; k++){
      outputfile << dt*k << "\t" << outVec_GS[k] << "\t" << outVec_ES[k] << "\t";
//       outputfile << Phase[k] << "\t" << Frequency[k];
      outputfile << endl;
    }
    outputfile.close();
    
//      outputfile.open("data/out_TS_full");
//      for(int k = k_tau; k < outVec_GS.size()-2; k++){
//        outputfile << dt*k << "\t" << outVec_GS[k] << "\t" << outVec_ES[k];
//  //       outputfile << "\t" << Phase[k] << "\t" << Frequency[k];
//        outputfile << "\t" << avfreq[k-k_tau] << "\t" << GSinv[k];
//  //       outputfile << "\t" << ESinv[k] << "\t" << outVecComplex_GS[k].real() << "\t" << outVecComplex_GS[k].imag();
//        outputfile << endl;
//      }
//      outputfile.close();
     
     
//      double GSinvAv = 0.0;
//      for(unsigned int k = 0; k < GSinv.size(); k++){
//        GSinvAv += GSinv[k];
//      }
//      cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
//      cout << eval_GS.average << "\t" << GSinvAv/GSinv.size() << endl;
//      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
       
    
//    //for Poincare Section
//    outputfile.open("data/out_TS_full_PS");
//    for(int k = k_tau; k < outVec_GS.size()-2; k++){
////       if(outVec_ES[k]/0.0835587 > 0.063 && outVec_ES[k]/0.0835587 < 0.067){
//      if(outVec_ES[k]/0.020171 > 0.063 && outVec_ES[k]/0.020171 < 0.067){
//        outputfile << dt*k << "\t" << outVec_GS[k] << "\t" << outVec_ES[k];
//  //       outputfile << "\t" << Phase[k] << "\t" << Frequency[k];
//        outputfile << "\t" << avfreq[k-k_tau] << "\t" << GSinv[k];
//  //       outputfile << "\t" << ESinv[k] << "\t" << outVecComplex_GS[k].real() << "\t" << outVecComplex_GS[k].imag();
//        outputfile << endl;
//      }
//    }
//    outputfile.close();
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-powerSpec" , false)){
      std::vector<double> powerSpec;
      std::vector<double> WindowedData;
      //GS
      sm::window_Hann(outVec_GS, WindowedData);
      katana::get_power_spec(WindowedData, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, "data/out_TS_powerSpec_GS_Hann");
      katana::sm_dump_power_spec(powerSpec, dt, 10, "data/out_TS_powerSpec_GS_Hann");
      //no window
      katana::get_power_spec(outVec_GS, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, "data/out_TS_powerSpec_GS");
      katana::sm_dump_power_spec(powerSpec, dt, 10, "data/out_TS_powerSpec_GS");
//       //ES
//       sm::window_Hann(outVec_ES, WindowedData);
//       katana::get_power_spec(WindowedData, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, "data/out_TS_powerSpec_ES_Hann");
//       //no window
//       katana::get_power_spec(outVec_ES, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50, "data/out_TS_powerSpec_ES");
      
//       std::vector<double> powerSpec;
//       katana::get_power_spec(outVec_GS, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50.0, "data/out_TS_powerSpec_GS");
//       katana::get_power_spec(outVec_ES, powerSpec);
//       katana::dump_power_spec(powerSpec, dt, powerSpec.size()/50.0, "data/out_TS_powerSpec_ES");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-opticalSpec" , false)){
      std::vector<double> opticalSpec;
      std::vector<std::complex<double>> WindowedData;
      sm::window_Hann(outVecComplex_GS, WindowedData);
      katana::get_optical_spec(WindowedData, opticalSpec);
//       katana::dump_optical_spec(opticalSpec, dt, opticalSpec.size()/10, "data/out_TS_opticalSpec_Hann");
      katana::sm_dump_optical_spec(opticalSpec, dt, 10, "data/out_TS_opticalSpec_Hann");
      //no window
      katana::get_optical_spec(outVecComplex_GS, opticalSpec);
//       katana::dump_optical_spec(opticalSpec, dt, opticalSpec.size()/20, "data/out_TS_opticalSpec");
      katana::sm_dump_optical_spec(opticalSpec, dt, 10, "data/out_TS_opticalSpec");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-AC" , false)){
      std::vector<double> AC;
      katana::get_autocorr_norm_out(outVecComplex_GS, AC);
      katana::print_correlation(AC, dt, "data/out_TS_AC_GS", AC.size(), 1, "ns");
      katana::get_autocorr_norm_out(outVecComplex_ES, AC);
      katana::print_correlation(AC, dt, "data/out_TS_AC_ES", AC.size(), 1, "ns");
    }
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
    
  }
  
  
  //K_GS sweep with CGS
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_KGS_wCGS" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.K_GS;
    str_sweep_parameter = "K_GS";
    PtrToSecPar = &p.C_GS;
    str_sweep_sec_parameter = "C_GS";
    
    str_sweep_updown = "up";
    sweep();

    if(bool_wdownsweep){
      sweep_start = sweep_end;
      sweep_step = -sweep_step;
      
      str_sweep_updown = "down";
      sweep();
    }
  }
  
  
  //K_GS sweep with CGS
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_wKGS" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J;
    str_sweep_parameter = "J";
    PtrToSecPar = &p.K_GS;
    str_sweep_sec_parameter = "K_GS";
    
    str_sweep_updown = "up";
    sweep();

    if(bool_wdownsweep){
      sweep_start = sweep_end;
      sweep_step = -sweep_step;
      
      str_sweep_updown = "down";
      sweep();
    }
  }
  


 
  return 0;
}


