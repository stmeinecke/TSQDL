struct parameters{
  //class to collect all parameters of the system.
  
  double J;		//pump current density
  double gGS;		//ground state gain
  double gES;		//excited state gain
  double kappaGS;	//ground state optical losses
  double kappaES;	//excited state optical lossen
  double beta;		//chance to spontaneously emit into lasing mode
  double ZQD;		//Number of quantum dots
  double NQD;		//area density of QDs
  double fact;		//fraction of active QDs
  double WGS;		//Einstein coefficient for spontaneous emission for the GS
  double WES;		//Einstein coefficient for spontaneous emission for the ES
  double RWloss;	//loss rate for the quantum wells
  double etaGS;		//electric field conversion factor for the GS
  double etaES;		//electric field conversion factor for the ES
  double GS_del_om_ES;		//GS frequency shift due to excited state population
  double GS_del_om_QW_e;	//GS frequency shift due to quantum well electron density
  double GS_del_om_QW_h;	//GS frequency shift due to quantum well hole density
  double ES_del_om_GS;		//GS frequency shift due to ground state population
  double ES_del_om_QW_e;	//GS frequency shift due to quantum well electron density
  double ES_del_om_QW_h;	//GS frequency shift due to quantum well hole density
  double gamma_GS;		//scales GS frequency shift due to carrier occupations
  double gamma_ES;		//scales ES frequency shift due to carrier occupations
  double img_offset;    //offset to the imaginare gain to center rotating frame
  
  double tau_GS;
  double tau_ES;
  double K_GS;
  double K_ES;
  double C_GS;
  double C_ES;
  
  double SqrtNoiseStr_GS;
  double SqrtNoiseStr_ES;
  std::complex<double> noise_GS;
  std::complex<double> noise_ES;
  
  double fakeSpE;
  bool bkick;
  double kick_time;
  double kick_strength;
  
  double SScale;

  
  parameters(){}
}; 
