class vars{ 
	 public:
		std::complex<double> EGS;	//complex ground state electric field amplitude in V * nm^-1
		std::complex<double> EES;	//complex excited state electric field amplitude in V * nm^-1
		
		double rhoGSe;		//electron ground state occupation probability
		double rhoGSh;		//hole ground state occupation probability
		double rhoESe;		//electron excited state occupation probability
		double rhoESh;		//hole excited state occupation probability
		
		double rhoGSei;		//inactive electron ground state occupation probability
		double rhoGShi;		//inactive hole ground state occupation probability
		double rhoESei;		//inactive electron excited state occupation probability
		double rhoEShi;		//inactive hole excited state occupation probability
		
		double we;			//quantum well electron area density in nm^-2
		double wh;			//qunatum well hole area density in nm^-2

		vars(){this->setTo(0.0);}

		void add(vars *v, vars *r){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i] + ((double*)v)[i];
			}
		} 

		void mult(vars *r, double m){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i]*m ;
			}
		}

		void setTo(double s){
			for(int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)this)[i] = s;
			}
		}

};

namespace IND{
  
  int EGS = 0;
  int EES = 2;
  int rhoGSe = 4;
  int rhoGSh = 5;
  int rhoESe = 6;
  int rhoESh = 7;
  int rhoGSei = 8;
  int rhoGShi = 9;
  int rhoESei = 10;
  int rhoEShi = 11;
  int we = 12;
  int wh = 13;
  
}
