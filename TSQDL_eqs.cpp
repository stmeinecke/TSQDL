
auto TSQDL = [](vars *x, vars_vec_wdX *Xhist, vars *d, parameters *p){
/// Scattering rates for T=300K, QD GS and ES
/// input we and wh in units of 10^11 cm^-2
/// Resulting rates are in units of 1/ns
    
  
    x->wh = x->we + (2.0*p->NQD*(p->fact*x->rhoGSe+(1.0-p->fact)*x->rhoGSei))
                  + (4.0*p->NQD*(p->fact*x->rhoESe+(1.0-p->fact)*x->rhoESei))
                  - (2.0*p->NQD*(p->fact*x->rhoGSh+(1.0-p->fact)*x->rhoGShi))
                  - (4.0*p->NQD*(p->fact*x->rhoESh+(1.0-p->fact)*x->rhoEShi));
  
    const double slowHCap = 1.0; //hole capture reduced to 5% required for ground state quenching. See Andres master thesis or understanding ground state quenching paper
    double S_GScap_e_in=0, S_EScap_e_in=0, S_GScap_h_in=0, S_EScap_h_in=0, S_GSESrel_e_in=0, S_GSESrel_h_in=0;
    double S_GScap_e_out=0, S_EScap_e_out=0, S_GScap_h_out=0, S_EScap_h_out=0, S_GSESrel_e_out=0, S_GSESrel_h_out=0;
    
    S_GScap_e_in= p->SScale * 18.5*x->we*x->we/(1.85+x->we);
    S_EScap_e_in= p->SScale * 48.3*x->we*x->we/(0.483+x->we);
    S_GScap_h_in= p->SScale * slowHCap * 10.5*x->wh*x->wh/(5.26+x->wh);
    S_EScap_h_in= p->SScale * slowHCap * 21.4*x->wh*x->wh/(1.79+x->wh);
    S_GSESrel_e_in= p->SScale * 1014*x->we/(1.39+x->we);
    S_GSESrel_h_in= p->SScale * 2272*x->wh/(2.27+x->wh);

    S_GScap_e_out=S_GScap_e_in*0.084133/(exp(0.215347*x->we)-1.0);
    S_EScap_e_out=S_EScap_e_in*0.581884/(exp(0.215347*x->we)-1.0);
    S_GScap_h_out=S_GScap_h_in*0.25828/(exp(0.0205776*x->wh)-1.0);
    S_EScap_h_out=S_EScap_h_in*0.559808/(exp(0.0205776*x->wh)-1.0);
    S_GSESrel_e_out=S_GSESrel_e_in*0.144587;
    S_GSESrel_h_out=S_GSESrel_h_in*0.461373;
    
  
    
  std::complex<double>i = std::complex<double>(0.0,1.0); // imaginary i
  const double saveSpE = 0E-6;
  
  //fields
  d->EGS = (p->gGS * (x->rhoGSe + x->rhoGSh - 1) - i * p->gamma_GS * ( p->GS_del_om_ES * (p->fact *(x->rhoESe + x->rhoESh) + (1-p->fact)*(x->rhoESei + x->rhoEShi)) + p->GS_del_om_QW_e * x->we + p->GS_del_om_QW_h * x->wh) - p->kappaGS) * x->EGS	//stimulated emission
#ifdef stoch_spon
	   + sqrt(x->rhoGSe * x->rhoGSh) * p->noise_GS	//spontaneous emisson
#endif  
#ifdef det_spon
	   + p->beta * p->ZQD * p->fact * p->WGS * p->etaGS * p->etaGS * x->rhoGSe * x->rhoGSh * x->EGS / (norm(x->EGS) + saveSpE)	//spontaneous emisson
#endif
     + i * p->img_offset * x->EGS	//offset to the imaginary gain to center rotating frame
// 	   + p->kappaGS * p->K_GS * exp(-i*p->C_GS) * Xhist->at_C(p->tau_GS, IND::EGS)
     + p->kappaGS * p->K_GS * exp(-i*p->C_GS) * Xhist->at_C_cubicHermite(p->tau_GS, IND::EGS)
     + p->fakeSpE;
#ifdef kick
  if(t < p->kick_time) d->EES += p->kick_strength;
#endif
	   
  
  d->EES = (p->gES * (x->rhoESe + x->rhoESh - 1) - i * p->gamma_ES * ( p->ES_del_om_GS * (p->fact *(x->rhoGSe + x->rhoGSh) + (1-p->fact)*(x->rhoGSei + x->rhoGShi)) + p->ES_del_om_QW_e * x->we + p->ES_del_om_QW_h * x->wh) - p->kappaES) * x->EES 	//stimulated emission
#ifdef stoch_spon
	   + sqrt(x->rhoESe * x->rhoESh) * p->noise_ES	//spontaneous emisson
#endif  
#ifdef det_spon
	   + p->beta * p->ZQD * p->fact * p->WES * p->etaES * p->etaES * x->rhoESe * x->rhoESh * x->EES / (norm(x->EES)+saveSpE) //spontaneous emission
#endif	   
//      + p->kappaES * p->K_ES * exp(-i*p->C_ES) * Xhist->at_C(p->tau_ES, IND::EES)
     + p->kappaES * p->K_ES * exp(-i*p->C_ES) * Xhist->at_C_cubicHermite(p->tau_ES, IND::EES)
     + p->fakeSpE;
#ifdef kick
  if(t < p->kick_time) d->EES += p->kick_strength;
#endif
  
  //active QD states
  d->rhoGSe = - (p->gGS / (p->ZQD * p->fact)) * (x->rhoGSe + x->rhoGSh - 1) * norm(x->EGS) / (p->etaGS * p->etaGS)
	      - p->WGS * x->rhoGSe * x->rhoGSh
	      + S_GScap_e_in * (1 - x->rhoGSe) - S_GScap_e_out * x->rhoGSe
	      + S_GSESrel_e_in * (1 - x->rhoGSe) * x->rhoESe - S_GSESrel_e_out * x->rhoGSe * (1 - x->rhoESe);

  d->rhoGSh = - (p->gGS / (p->ZQD * p->fact)) * (x->rhoGSe + x->rhoGSh - 1) * norm(x->EGS) / (p->etaGS * p->etaGS)
	      - p->WGS * x->rhoGSe * x->rhoGSh
	      + S_GScap_h_in * (1 - x->rhoGSh) - S_GScap_h_out * x->rhoGSh
	      + S_GSESrel_h_in * (1 - x->rhoGSh) * x->rhoESh - S_GSESrel_h_out * x->rhoGSh * (1 - x->rhoESh);
	      
  d->rhoESe = - (p->gES / (2 * p->ZQD * p->fact)) * (x->rhoESe + x->rhoESh - 1) * norm(x->EES) / (p->etaES * p->etaES)
	      - p->WES * x->rhoESe * x->rhoESh
	      + S_EScap_e_in * (1 - x->rhoESe) - S_EScap_e_out * x->rhoESe
	      - 0.5 * (S_GSESrel_e_in * (1 - x->rhoGSe) * x->rhoESe - S_GSESrel_e_out * x->rhoGSe * (1 - x->rhoESe));
	      
  d->rhoESh = - (p->gES / (2 * p->ZQD * p->fact)) * (x->rhoESe + x->rhoESh - 1) * norm(x->EES) / (p->etaES * p->etaES)
	      - p->WES * x->rhoESe * x->rhoESh
	      + S_EScap_h_in * (1 - x->rhoESh) - S_EScap_h_out * x->rhoESh
	      - 0.5 * (S_GSESrel_h_in * (1 - x->rhoGSh) * x->rhoESh - S_GSESrel_h_out * x->rhoGSh * (1 - x->rhoESh));
	      
  //inactive QDs      
  d->rhoGSei = - p->WGS * x->rhoGSei * x->rhoGShi
	      + S_GScap_e_in * (1 - x->rhoGSei) - S_GScap_e_out * x->rhoGSei
	      + S_GSESrel_e_in * (1 - x->rhoGSei) * x->rhoESei - S_GSESrel_e_out * x->rhoGSei * (1 - x->rhoESei);

  d->rhoGShi = - p->WGS * x->rhoGSei * x->rhoGShi
	      + S_GScap_h_in * (1 - x->rhoGShi) - S_GScap_h_out * x->rhoGShi
	      + S_GSESrel_h_in * (1 - x->rhoGShi) * x->rhoEShi - S_GSESrel_h_out * x->rhoGShi* (1 - x->rhoEShi);      
	      
  d->rhoESei = - p->WES * x->rhoESei * x->rhoEShi
	      + S_EScap_e_in * (1 - x->rhoESei) - S_EScap_e_out * x->rhoESei
	      - 0.5 * (S_GSESrel_e_in * (1 - x->rhoGSei) * x->rhoESei - S_GSESrel_e_out * x->rhoGSei * (1 - x->rhoESei));
	      
  d->rhoEShi = - p->WES * x->rhoESei * x->rhoEShi
	      + S_EScap_h_in * (1 - x->rhoEShi) - S_EScap_h_out * x->rhoEShi
	      - 0.5 * (S_GSESrel_h_in * (1 - x->rhoGShi) * x->rhoEShi - S_GSESrel_h_out * x->rhoGShi * (1 - x->rhoEShi));
	           
  //carrier area densities
  d->we = p->J - p->RWloss * x->we * x->wh
	  - 2 * p->NQD * p->fact * (S_GScap_e_in * (1 - x->rhoGSe) - S_GScap_e_out * x->rhoGSe)
	  - 4 * p->NQD * p->fact * (S_EScap_e_in * (1 - x->rhoESe) - S_EScap_e_out * x->rhoESe)
	  - 2 * p->NQD * (1 - p->fact) * (S_GScap_e_in * (1 - x->rhoGSei) - S_GScap_e_out * x->rhoGSei)
	  - 4 * p->NQD * (1 - p->fact) * (S_EScap_e_in * (1 - x->rhoESei) - S_EScap_e_out * x->rhoESei);
	  
//   d->wh = p->J - p->RWloss * x->we * x->wh
// 	  - 2 * p->NQD * p->fact * (S_GScap_h_in * (1 - x->rhoGSh) - S_GScap_h_out * x->rhoGSh)
// 	  - 4 * p->NQD * p->fact * (S_EScap_h_in * (1 - x->rhoESh) - S_EScap_h_out * x->rhoESh)
// 	  - 2 * p->NQD * (1 - p->fact) * (S_GScap_h_in * (1 - x->rhoGShi) - S_GScap_h_out * x->rhoGShi)
// 	  - 4 * p->NQD * (1 - p->fact) * (S_EScap_h_in * (1 - x->rhoEShi) - S_EScap_h_out * x->rhoEShi);
	  	  
};
