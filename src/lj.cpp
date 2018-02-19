/**
//  lj.cpp
//  md_modeling
//
//  Created by Swabir Silayi
//  Copyright © 2018 Swabir Silayi. All rights reserved.
*/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

int N = 256; 		         /*number of atoms*/
int nSteps = 20000;	     /*number of time-steps*/
double dt = 0.001;           /*time step*/
double T ;                   /*temperature*/
double rho ;                 /*number density*/


double **r;                 /*positions*/
double **v;                 /*velocities*/
double **f;                 /*forces*/

double a;                   /*lattice constant/ unit cell size*/
double L;                   /*box size*/
double rCutOff;             /*cut off for interactions*/
double vol;                 /*box volume*/


double U;                   /*potential*/
double K;                   /*kinetic*/
double E;                   /*total energ*/
double vir;                 /*virial*/



/*variables for the g(r) function*/
double * gr;                /*radial distribution bins*/
int nbins;                  /*number of bins*/
double dbin = 0.025;          /*size of bins*/


/*variables for the velocity correlation function*/
double sum0;
double sum;
double **v_old;

/*initialize positions in fcc lattice structure with random velocities*/
void initialize( )
{
    r = new double* [N];
    v = new double* [N];
    f = new double* [N];
    
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        f[i] = new double [3];
    }
    
    
   
    //L from N and rho
    L = pow(N / rho, 1.0/3);
    
    //find M for N atoms in fcc lattice
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    
    //lattice constant and cut off
    a = L / M;    
    rCutOff = 0.49*L;
    
    //volume
    vol = L*L*L;

   
    // positions in fcc unit cell
    double dx[4] = {0.0, 0.5, 0.5, 0.0};
    double dy[4] = {0.0, 0.5, 0.0, 0.5};
    double dz[4] = {0.0, 0.0, 0.5, 0.5};
    
    
    int n = 0;
    
    for (int x = 0; x < M; x++)
        for (int y = 0; y < M; y++)
            for (int z = 0; z < M; z++)
                for (int p = 0; p < 4; p++){
                    if (n < N){
                        r[n][0] = (x + dx[p]) * a;
                        r[n][1] = (y + dy[p]) * a;
                        r[n][2] = (z + dz[p]) * a;
                        ++n;
                    }
                }
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            v[i][j] = rand()-0.5;
    
    //center of mass positions and velocities
    double rCM[3] = {0, 0, 0};
    double vCM[3] = {0, 0, 0};
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
        {
            rCM[j] += r[i][j];
            vCM[j] += v[i][j];
        }
    
    for (int i = 0; i < 3; i++)
    {   rCM[i] /= double(N);
        vCM[i] /= double(N);
    }
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
        {
            r[i][j] -= rCM[j];
            v[i][j] -= vCM[j];
            
        }
    
    
    //rescale velocities to temperature T
    double vSqdSum = 0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            vSqdSum += v[i][j] * v[i][j];
    
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            v[i][j] *= lambda;
    
    /*radial distribution function variablee initialization*/
    //bin count for radial distribution
    nbins =  (int)( rCutOff / dbin) + 1 ;
    gr = new double [nbins];
    
    /*velocity autocorrelation funtion variable initialization*/
    v_old = new double* [N];
    
    for (int i = 0; i < N; i++) {
        v_old[i] = new double [3];
    }
     sum0 = 0.0;
     for (int i = 0; i < N; i++){
     
        v_old [i][0] = v[i][0];
        v_old [i][1] = v[i][2];
        v_old [i][2] = v[i][2];
        
        sum0 = sum0 + v_old[i][0]* v_old[i][0] + v_old[i][1]* v_old[i][1] + v_old[i][2]* v_old[i][2];
        
    }
    
    
    
    return;
    
}

/*velocity rescaling to target temperature*/
void rescaleVelocities() {
    
    double vSqdSum = 0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            vSqdSum += v[i][j] * v[i][j];
    
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            v[i][j] *= lambda;
  
}


/*force calculation*/
void force ( ){
    
    U = 0.0;
    vir = 0.0;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 3; j++)
            f[i][j] = 0.0;
    
    double dr[3] = {0.0}; double rSqd = 0.0;
    
    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            
            dr[0] = r[i][0] - r[j][0];
            dr[1] = r[i][1] - r[j][1];
            dr[2] = r[i][2] - r[j][2];
            
            /*minimum image convention*/
            for (int k = 0; k < 3; k++){
                if (dr[k] > 0.5*L)
                    dr[k] = dr[k]- L;
                
                if (dr[k] < -0.5*L)
                    dr[k] = dr[k]+L;
            }
            
            
            rSqd = dr[0]*dr[0]+ dr[1]*dr[1]+ dr[2]*dr[2];
            
            if (rSqd <= rCutOff*rCutOff)
            {
                
                
                double ff = 24.0* (2.0/ pow(rSqd, 7) - 1.0/ pow(rSqd, 4));
                
                double fc = 24.0* (2.0/ pow((rCutOff*rCutOff), 7) - 1.0/ pow((rCutOff*rCutOff), 4));
                
                ff = ff - fc;
                
                f[i][0] = f[i][0] + (ff* dr[0] );
                f[j][0] = f[j][0] - (ff* dr[0] );
                
                f[i][1] = f[i][1] + (ff* dr[1]) ;
                f[j][1] = f[j][1] - (ff* dr[1]) ;
                
                f[i][2] = f[i][2] + (ff* dr[2] );
                f[j][2] = f[j][2] - (ff* dr[2]) ;
                
                double U_cut =  4.0*(1.0/ pow((rCutOff*rCutOff), 6) - 1.0/ pow((rCutOff*rCutOff), 3));
                double Uf_cut = sqrt(rSqd) - rCutOff;
                
                Uf_cut = Uf_cut*fc;
                U_cut = U_cut + Uf_cut;
                
                U = U + 4.0*(1.0/ pow(rSqd, 6) - 1.0/ pow(rSqd, 3)) - U_cut;
                
                
                vir += ff*rSqd;
            }
        }//end for j
    }//end for i
    
    return;
}


/*record positions and velocities to in .xyz format*/
void record(string filename){
    
    ofstream file; file.open(filename.c_str());
    
    file << N <<"\n";
    file << L <<"\n";
    for (int i = 0; i < N; i++)
    {
        file << "Ar" <<"\t";
        for (int j = 0; j < 3; j++)
            file << r[i][j] << setw(12) <<"\t";
        for (int j = 0; j < 3; j++)
            file << v[i][j] << setw(12) <<"\t";
        file << "\n";
    }
    
    file.close();
}

/*quantity running average and variance calculation*/
void ave_var(int cc, double *UAvg, double *Var, double UTOTAL)
{
    
    UAvg[cc] = cc == 0 ? UTOTAL : UAvg[cc - 1] + ((UTOTAL - UAvg[cc - 1]) / cc);
    
    Var[cc] = cc == 0 ? 0 : Var[cc - 1] * cc + (UTOTAL - UAvg[cc - 1]) * (UTOTAL - UAvg[cc]);
    Var[cc] /= (cc + 1);
}


/*radial distribution function*/
void update_gr ( ) {
    
    double dr[3] = {0.0}; double rSqd = 0.0;
    
    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            
            dr[0] = r[i][0] - r[j][0];
            dr[1] = r[i][1] - r[j][1];
            dr[2] = r[i][2] - r[j][2];
            
            
            for (int k = 0; k < 3; k++){
                if (dr[k] > 0.5*L)
                    dr[k] = dr[k]- L;
                
                if (dr[k] < -0.5*L)
                    dr[k] = dr[k]+L;
            }
            
            rSqd = dr[0]*dr[0]+ dr[1]*dr[1]+ dr[2]*dr[2];
            
            if (rSqd <= rCutOff*rCutOff)
            {
                int bin=(int)(sqrt(rSqd)/dbin);
                gr[bin]+=2;
            }
        }
    }
    
    return;
}


/*translational order parameter*/
double order_parameter( ){
    
    double lmbd = 0.0;
    
    double lx = 0.0, ly = 0.0, lz = 0.0;
    
    for (int i = 0; i < N; i++)
    {
        
        lx = lx + cos ((4.0*M_PI*r[i][0])/a);
        ly = ly + cos ((4.0*M_PI*r[i][1])/a);
        lz = lz + cos ((4.0*M_PI*r[i][2])/a);
    }
    
    lx = lx/double (N);
    ly = ly/double (N);
    lz = lz/double (N);
    
    lmbd = (1.0/3.0)*(lx + ly + lz);
    
    return lmbd;
    
}

/*velocity autocorrelation function*/
double vacf(){
    sum = 0.0;
    for (int i = 0; i < N; i++){
    
      sum = sum +   v[i][0]*v_old[i][0] +  v[i][1]*v_old[i][1] +  v[i][2]*v_old[i][2];
        
    }
    
    double C = sum/sum0;
    
    return C;

}

int main( ){
    
    
    
    double ti[11] = {2.0,1.5,1.0,0.8,0.75,0.7,0.65,0.6,0.5,0.3,0.1};
    double ri[2] = {1.0, 0.7};
    
    for (int rn = 0; rn < 2; rn ++){
        
        rho  = ri[rn];
        
        for (int tn = 0; tn < 11; tn++){
            
            T = ti[tn];
            
            if (T == 1.5 || T == 0.1){ 
                
                /*initialize positions and velocities and calculate the initial forces*/ 
                initialize( ); force( );
                cout << U << "\n";
                
                cout << "running for rho: "<< rho << ", T : "<< T << "\n"; 
                
                
                
                
                /*make space for average values*/
                int nave = int(nSteps);
                double * U_ave = new double[nave];
                double * K_ave = new double [nave];
                double * U_var = new double [nave];
                double * K_var = new double [nave];
                
                double * E_ave = new double [nave];
                double * E_var = new double [nave];
                
                double * t_ave = new double [nave];
                double * t_var = new double [nave];
                
                double * p_ave = new double [nave];
                double * p_var = new double [nave];
                
                
                /*open files and write column headers to file*/
                std::ostringstream fn;
                fn << rho <<"/init_"<<T<<".xyz";
                std::string fn1 = fn.str();
                record (fn1.c_str());
                
                std::ostringstream in;
                in <<rho<< "/info_"<<T<<".txt";
                std::string in1 = in.str();
                ofstream file; file.open (in1.c_str());
                
                std::ostringstream gn;
                gn << rho<<"/rdf_"<<T<<".txt";
                std::string gn1 = gn.str();
                ofstream gfile; gfile.open (gn1.c_str());
                
                std::ostringstream vn;
                vn << rho<<"/vacf_"<<T<<".txt";
                std::string vn1 = vn.str();
                ofstream vfile; vfile.open (vn1.c_str());
                
                
                
                std::ostringstream lmf;
                lmf << rho<< "/lmf_"<<T<<".txt";
                std::string ln1 = lmf.str();
                ofstream lfile; lfile.open (ln1.c_str());
                
                
                
                file << "timestep "<<setw(12) << "\t" 
                << "T   "<<setw(12) << "\t" << "T_{ave}  "<<setw(12) << "\t" << "T_{var}"<<setw(12) << "\t" 
                << "K/N "<<setw(12) << "\t" << "K_{ave}/N"<<setw(12) << "\t" << "K_{var}"<<setw(12) << "\t" 
                << "U/N "<<setw(12) << "\t" << "U_{ave}/N"<<setw(12) << "\t" << "U_{var}"<<setw(12) << "\t" 
                << "E/N "<<setw(12) << "\t" << "E_{ave}/N"<<setw(12) << "\t" << "E_{var}"<<setw(12) << "\t" 
                << "P   "<<setw(12) << "\t" << "P_{ave}  "<<setw(12) << "\t" << "P_{var}" <<"\n";
                
                gfile << "Distance(r)"<<setw(12) << "\t" << "(g[r])"
                << "\t" << "(gCs[r])" <<setw(12)
                << "\t" << "(gK[r])"  <<setw(12)
                << "\t" << "(gCsK[k])"  <<"\n";
                
                
                
                lfile <<  "time-step" << "\t" << "lambda" <<"\n";
                vfile <<  "time-step" << "\t" << "Vacf(t)" <<"\n";
                vfile <<  0 << "\t" << vacf() <<"\n";
                
                
                
                /*loop over time*/
                
                int cc = 0; int ngr = 0; int vcnt = 0;
                time_t si = time (NULL);
                
                
                for (int n = 0; n < nSteps; n++)
                {
                    
                    //if (n == 0 || n%200 == 0){  
                        
                      vfile <<  n+1 << "\t" << vacf() <<"\n";
                        /*save configuration*/
                   /*     std::ostringstream vin;
                        vin << "vfiles/"<< rho <<"/"<< T << "/v_"<<vcnt<<".xyz";
                        std::string vin1 = vin.str();
                        record (vin1.c_str());
                        vcnt++;
                   // }
                    */
                    /*start velocity verlet integration scheme*/
                    /*loop over atoms to update positions*/
                    for (int i = 0; i < N; i++){
                        for (int k = 0; k < 3; k++){
                            
                            r[i][k] = r[i][k] + v[i][k]*dt + 0.5* f[i][k]* dt*dt;
                            
                            
                            /* apply periodic boundary conditions*/
                            if (r[i][k] > L)
                                r[i][k] = r[i][k] - L;
                            
                            if (r[i][k] < 0.0)
                                r[i][k] = r[i][k] + L;
                            
                            v[i][k] = v[i][k] + 0.5 * f[i][k]*dt;
                            
                            
                        }
                        
                    }/*end loop over atoms*/
                    
                    /*update force*/
                    force ( );
                    
                    /*update kinetic energy and velocites*/
                    K = 0.0;
                    
                    /*loop over atoms*/
                    for (int i = 0; i < N; i++){
                        for (int j = 0; j < 3; j++){
                            
                            v[i][j] = v[i][j] + 0.5* f[i][j]*dt;
                            
                            K = K + 0.5 * v[i][j] * v[i][j];
                        }
                    }/*end loop over atoms*/
                    
                    
                    /*total energy*/
                    E = K+U;
                    
                    /*end velocity verlet integration scheme*/
                    
                    
                    /*production steps after completing the equlibration steps*/
                    if (n > 0.5*nSteps){
                        
                        
                        /*to calculate system temperature and velocity autocorrelation*/
                        double vSqdSum = 0;
                        
                        for (int i = 0; i < N; i++)
                            for (int j = 0; j < 3; j++){
                                vSqdSum += v[i][j] * v[i][j];
                                
                                
                                
                            }
                            
                            
                            double temp = ( vSqdSum/ (3*(N-1)) );
                        
                        /*system pressure from Kinetic energy and virial*/
                        double P = rho*temp + vir*(rho/(3*(N-1)));
                        
                        
                        //calculate averages
                        ave_var(cc, U_ave, U_var, U);
                        ave_var(cc, K_ave, K_var, K);
                        ave_var(cc, E_ave, E_var, E);
                        ave_var(cc, t_ave, t_var, temp);
                        ave_var(cc, p_ave, p_var, P);
                        
                        
                        /*record averages every 200 time steps*/
                        if (cc %200 == 0)
                        {
                            update_gr ( ); ngr++;
                            
                            file << n << "\t" << temp << "\t" << t_ave[cc] << "\t" << t_var[cc] << "\t"
                            << K / double (N) << "\t" << K_ave[cc] / double (N) << "\t" << K_var[cc] << "\t"
                            << U / double (N) << "\t" << U_ave[cc] / double (N) << "\t" << U_var[cc] << "\t"
                            << E / double (N) << "\t" << E_ave[cc] / double (N) << "\t" << E_var[cc] << "\t"
                            << P << "\t" << p_ave[cc] << "\t" << p_var[cc] << "\n";
                            
                        }
                        cc = cc+1;
                        
                        
                        
                    }/*end production loop*/
                    
                    
                    /*calculate the positional order parameter and save to file*/
                    lfile << n << "\t" << order_parameter() << "\n";
                    
                    
                    /*rescale velocities*/
                    /*    if (n % 100 == 0)
                     *                    rescaleVelocities();
                     */
                    
                }//end loop over time
                
                file.close();
                
                
                time_t sf = time(NULL);
                
                
                /*save final configuration*/
                std::ostringstream fin;
                fin << rho << "/fin_"<<T<<".xyz";
                std::string fin1 = fin.str();
                record (fin1.c_str());
                
                
                
                /* Normalize radial distribution g(r) and save  to file*/
                for (int i=0; i< nbins;i++) {
                    double rr = dbin*(i+0.5);
                    double vb=((i+1)*(i+1)*(i+1)-i*i*i)*dbin*dbin*dbin;
                    double nid=(4./3.)*M_PI*vb*rho;
                    
                    gfile << i*dbin<<"\t" <<(double)(gr[i])/(ngr*N*nid) <<"\n";
                }
                
                
                /* compute the long range corrections */
                double r3 = 1.0/(rCutOff*rCutOff*rCutOff);
                double ucor = -8*M_PI*rho*r3/3.0;
                
                double vSqdSum = 0;    
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < 3; j++)
                        vSqdSum += v[i][j] * v[i][j];
                    
                    double temp = ( vSqdSum/(3 * (N-1)) );
                
                double pcor = -16.0/3.0*M_PI*rho*r3/(3.0*temp);
                
                
                /*record final averge quantities*/
                std::ostringstream trn;
                trn << rho<<"/tracking.txt";
                std::string trn1 = trn.str();
                ofstream tfile; tfile.open (trn1.c_str(), ios_base::app);
                
                
                
                tfile << T << "\t" << K_ave[cc-1] / double (N) << "\t" << U_ave[cc-1] / double (N) << "\t" 
                << E_ave[cc-1]/ double (N) << "\t" << p_ave[cc-1]<< "\t" << t_ave[cc-1]
                << "\t U_corr : " << ucor <<"\t P_corr : " << pcor <<"\t runtime:" << sf-si << "\n";
                
                
                
                file.close();
                gfile.close();
                vfile.close();
                lfile.close();
                tfile.close();
                
                
            }
        }/*end for loop over T*/
        
    }/*end for loop over rho*/
    
    
    return 0;
    
    
}

/*EOF*/

