



class evaluateTS {
public:
  double average, greatestMax, smallestMin, period;
  double dt;
  std::vector<double> maxima;
  std::vector<double> minima;
  std::vector<double> uniqueMax;
  std::vector<double> uniqueMin;
  std::vector<double>* datavec;

  evaluateTS();
  evaluateTS(std::vector<double>* datavec, double dt);
  void newTS(std::vector<double>* datavec, double dt);  
  int numberOfMax(void);
  int numberOfMin(void);
  bool checkTransient(void);
  double opticalFreq(std::vector<std::complex<double>> &ComplexE, std::vector<double> &Phase, std::vector<double> &Frequency, double dt);
  double loopsRatio(void);
  
  double extremaTol = 1E-16;
  double varTol = 0.5E-3;
  double avTol = 0.05; //tolarance for detecting a drifting average for transient detection
  double doubleCountTol = 0.05;
  double max_min_dscr = 1E-5; //used as a cut-off for detecting limit cycles/extrema by relativ difference between greatest max and smallest min
  double av_thresh = 1E-6;

  
private:
  int size;
  void basics();
  double t_i, t_f;
  double interpolateQuadMax(const double X1, const double f1, const double X2, const double f2, const double X3, const double f3);
};
 
evaluateTS::evaluateTS(){};

evaluateTS::evaluateTS(std::vector<double>* datavec, double dt){
  this->size = datavec->size();
  this->datavec = datavec;
  this->dt = dt;
  basics();
}

void evaluateTS::newTS(std::vector<double>* datavec, double dt){
  this->size = datavec->size();
  this->datavec = datavec;
  this->dt = dt;
  basics();
}

void evaluateTS::basics(){
  maxima.resize(0);
  minima.resize(0);
  uniqueMax.resize(0);
  uniqueMin.resize(0);
  average = 0;
  greatestMax = datavec->at(0);
  smallestMin = datavec->at(0);
  period = 0;
  double sum = 0;;
  smallestMin = datavec->at(0);

  for(int k=0; k < size-2; k++){
    //calc sum for average
    sum += datavec->at(k);
    //find greates maximum
    if(greatestMax < datavec->at(k)){
      greatestMax = datavec->at(k);
    }
    //find smallest mininum
    if(smallestMin > datavec->at(k)){
      smallestMin = datavec->at(k);
    }
    //find maxima
    if(datavec->at(k) < (1.0-extremaTol) * datavec->at(k+1) && (1.0-extremaTol)*datavec->at(k+1) > datavec->at(k+2)){
      if(maxima.size() == 0) t_i = this->dt*(k);
      else t_f = this->dt*(k);
      //interpolate true maximum from with next neighbors 
      maxima.push_back(interpolateQuadMax(this->dt*(k),datavec->at(k),this->dt*(k+1),datavec->at(k+1),this->dt*(k+2),datavec->at(k+2)));
    }
    //find minima
    if(datavec->at(k) > (1.0+extremaTol) * datavec->at(k+1) && (1.0+extremaTol)*datavec->at(k+1) < datavec->at(k+2)){
      //interpolate true maximum from with next neighbors 
      minima.push_back(interpolateQuadMax(this->dt*(k),datavec->at(k),this->dt*(k+1),datavec->at(k+1),this->dt*(k+2),datavec->at(k+2)));
    }
  }
  average = sum / (size-2);
  if(maxima.size() > 1 && greatestMax/smallestMin > (1.0 + max_min_dscr) && average > av_thresh) period = (t_f - t_i)/maxima.size();
  else period = -1;
}

int evaluateTS::numberOfMax(){
//   if(average > 0.0 && (greatestMax/smallestMin < (1.0 + max_min_dscr) || fabs(average) < av_thresh) ) return 0;
//   else if(average < 0.0 && (smallestMin/greatestMax < (1.0 + max_min_dscr) || fabs(average) < av_thresh) ) return 0;
  if(average > 0.0 && ((greatestMax-smallestMin)/fabs(average) <  + max_min_dscr || fabs(average) < av_thresh) ) return 0;
  else if(average < 0.0 && ((greatestMax-smallestMin)/fabs(average) <  + max_min_dscr || fabs(average) < av_thresh) ) return 0;
  else{
    int max_counter = 0;
    bool doublecount;
    for(int k = 0; k < maxima.size(); k++){
      doublecount = false;
      for(int l = 0; l < max_counter; l++){
        if(fabs(maxima[k]-average) > (1.0-doubleCountTol) * fabs(uniqueMax[l]-average) && fabs(maxima[k]-average) < (1.0+doubleCountTol) * fabs(uniqueMax[l]-average)){
          doublecount = true;
        }
      }
      if(doublecount == false){
        uniqueMax.push_back(maxima[k]);
        max_counter += 1;
      }
    }
    return max_counter;
  }
}

int evaluateTS::numberOfMin(void){
//   if(average > 0.0 && (greatestMax/smallestMin < (1.0 + max_min_dscr) || fabs(average) < av_thresh) ) return 0;
//   else if(average < 0.0 && (smallestMin/greatestMax < (1.0 + max_min_dscr) || fabs(average) < av_thresh) ) return 0;
  if(average > 0.0 && ((greatestMax-smallestMin)/fabs(average) <  + max_min_dscr || fabs(average) < av_thresh) ) return 0;
  else if(average < 0.0 && ((greatestMax-smallestMin)/fabs(average) <  + max_min_dscr || fabs(average) < av_thresh) ) return 0;
  else{
    int min_counter = 0;
    bool doublecount;
    for(int k = 0; k < minima.size(); k++){
      doublecount = false;
      for(int l = 0; l < min_counter; l++){
        if(fabs(minima[k]-average) > (1.0-doubleCountTol) * fabs(uniqueMin[l]-average) && fabs(minima[k]-average) < (1.0+doubleCountTol) * fabs(uniqueMin[l]-average)){
          doublecount = true;
        }
      }
      if(doublecount == false){
        uniqueMin.push_back(minima[k]);
        min_counter += 1;
      }
    }
    return min_counter;
  }
}

double evaluateTS::interpolateQuadMax(const double X1, const double f1, const double X2, const double f2, const double X3, const double f3){
  //interpolates true extremum where f2 is a local extremum at X2 in the datavec with neighboring points f1 and f3
  const double x1 = 0.0;
  const double x2 = X2 - X1;
  const double x3 = X3 - X1;
  
  const double k1 = f1/((x1-x2)*(x1-x3));
  const double k2 = f2/((x2-x1)*(x2-x3));
  const double k3 = f3/((x3-x1)*(x3-x2));
  
  const double a = k1 + k2 + k3;
  const double b = -(k1*(x2+x3) + k2*(x1+x3) + k3*(x1+x2));
  const double c = k1*x2*x3 + k2*x1*x3 + k3*x1*x2;
  
  const double xmax = -0.5*b/a+X1;
  const double fmax = 0.25*b*b/a - 0.5*b*b/a + c;
  
  return fmax;
}

bool evaluateTS::checkTransient(void){
  //checks a datavec for either decay or growths by comparing the variance of the extrema in three consecutive sections of the datavec.
  if(greatestMax/smallestMin < (1.0 + max_min_dscr) || average < av_thresh){
    return false;
  }
  else{
    //calculate averages for two sections of the datavec
    double average1, average2;
    average1 = average2 = 0;
    for(int k = 0; k < size/2; k++){
      average1 += datavec->at(k);
    }
    average1 = 2*average1/size;
    for(int k = size/2; k < size - 1; k++){
      average2 += datavec->at(k);
    }
    average2 = 2*average2/size;
    //check for drifting overall average
    if(average1 < average2*(1.0-avTol) || average1 > average2*(1.0+avTol)) return true;
      
    if(maxima.size() > 30 && minima.size() > 30){
      //calculate extrema averages from the extrema for the three sections
      double av1, av2, av3;
      av1 = av2 = av3 = 0.0;
      int n_extrema = std::min((int)maxima.size(), (int)minima.size());
      int i_ext = 0;
      int N = 0;
      
      while(i_ext<n_extrema/3){
	av1 += maxima[i_ext] + minima[i_ext];
	i_ext += 1;
	N += 1;
      }
      av1 = av1/(2*N);
      N = 0;
      while(i_ext<2*n_extrema/3){
	av2 += maxima[i_ext] + minima[i_ext];
	i_ext += 1;
	N += 1;
      }
      av2 = av2/(N*2);
      N = 0;
      while(i_ext<n_extrema){
	av3 += maxima[i_ext] + minima[i_ext];
	i_ext += 1;
	N += 1;
      }
      av3 = av3/(N*2);
      
      //check for drifting extrama average
      if((av1 > (1.0+varTol)*av2 && av2 > (1.0+varTol)*av3) || (av1 < (1.0-varTol)*av2 && av2 < (1.0-varTol)*av3)) return true;
      else{
	//calculate the variances for the three sections using the previouly calculated averages
	double var1, var2, var3;
	var1 = var2 = var3 = 0.0;
	i_ext = 0;
	N = 0;
	while(i_ext < n_extrema/3){
	  var1 += (av1 - maxima[i_ext])*(av1 - maxima[i_ext]) + (av1-minima[i_ext])*(av1-minima[i_ext]);
	  i_ext += 1;
	  N += 1;
	}
	var1 = var1 / (2*N);
	N = 0;
	while(i_ext < 2*n_extrema/3){
	  var2 += (av2 - maxima[i_ext])*(av2 - maxima[i_ext]) + (av2-minima[i_ext])*(av2-minima[i_ext]);
	  i_ext += 1;
	  N += 1;
	}
	var2 = var2 / (2*N);
	N = 0;
	while(i_ext < n_extrema){
	  var3 += (av3 - maxima[i_ext])*(av3 - maxima[i_ext]) + (av3-minima[i_ext])*(av3-minima[i_ext]);
	  i_ext += 1;
	  N += 1;
	}
	var3 = var3 / (2*N);
	
// 	std::cout << var1 << "\t" << var2 << "\t" << var3 << std::endl;
// 	std::ofstream f;
// 	f.open("output/checkTransient");
// 	f << "# " << maxima.size() << "\t " << minima.size() << std::endl;
// 	int n=0;
// 	while(n < n_extrema/3){
// 	  f << n << "\t" << maxima[n] << "\t" << minima[n] << "\t" << av1 << "\t" << var1 << std::endl;
// 	  n += 1;
// 	}
// 	while(n < 2*n_extrema/3){
// 	  f << n << "\t" << maxima[n] << "\t" << minima[n] << "\t" << av2 << "\t" << var2 << std::endl;
// 	  n += 1;
// 	}
// 	while(n < n_extrema){
// 	  f << n << "\t" << maxima[n] << "\t" << minima[n] << "\t" << av3 << "\t" << var3 << std::endl;
// 	  n += 1;
// 	}
// 	f.close();
// 	
// 	std::cout << maxima.size() << " " << minima.size() << std::endl;
	if((var1 > (1.0+varTol)*var2 && var2 > (1.0+varTol)*var3) || (var1 < (1.0-varTol)*var2 && var2 < (1.0-varTol)*var3)) return true;
	else return false;
      }
    }
    else return false;
  }
}

// double evaluateTS::opticalFreq(){
//   //use with phase vector
//   int npi = 0;
//   for(int k = 0; k < size - 1; k++){
//     datavec->at(k) = datavec->at(k) + npi * sm::PI;
//     if (datavec->at(k) > sm::PI + datavec->at(k+1) + npi * sm::PI){
//       npi += 2;
//     }
//     else if (datavec->at(k) < -sm::PI + datavec->at(k+1) + npi * sm::PI){
//       npi -= 2;
//     }
//   }
//   double OptFreqAv = 0;
//   for(int k = 0; k < size-2;k++){
//     OptFreqAv = OptFreqAv + datavec->at(k+1) - datavec->at(k);
//   }
//   return OptFreqAv / (tvec->at(tvec->size()-2) - tvec->at(0));
// }


double evaluateTS::opticalFreq(std::vector<std::complex<double>> &ComplexE, std::vector<double> &Phase, std::vector<double> &Frequency, double dt){
  Phase.resize(ComplexE.size());
  Frequency.resize(ComplexE.size());
  int npi = 0;
  
  for(int k = 0; k < ComplexE.size() - 1; k++){
    Phase[k] = arg(ComplexE[k]) + npi * sm::PI;
    if (Phase[k] > sm::PI + arg(ComplexE[k+1]) + npi * sm::PI){
      npi += 2;
    }
    else if (Phase[k] < -sm::PI + arg(ComplexE[k+1]) + npi * sm::PI){
      npi -= 2;
    }
  }
  double OptFreqAv = 0;
  for(int k = 0; k < size-2;k++){
    Frequency[k] = (Phase[k+1] - Phase[k])/dt;
    OptFreqAv += Frequency[k];
  }
  return OptFreqAv / (ComplexE.size()-2);
}

double evaluateTS::loopsRatio(void){
  double lowerLoopMax = smallestMin + 0.75*(greatestMax-smallestMin);
  double upperLoopMin = smallestMin + 0.25*(greatestMax-smallestMin);
  
  int counterLLmaxima = 0;
  int counterULmininma = 0;
  for(unsigned int k = 0; k < maxima.size(); k++){
    if(maxima[k] < lowerLoopMax) counterLLmaxima ++;
  }
  for(unsigned int k = 0; k < minima.size(); k++){
    if(minima[k] > upperLoopMin) counterULmininma ++;
  }
  
  return (double)counterULmininma/(double)(counterLLmaxima+counterULmininma);
}
