/* 
This programm allows to determine the osmotic angular velocity \Omega^u_\vartheta
as a function of \vartheta. The code uses the numerical algorithm 
given in the Appendix of the thesis "Particle Spin described by 
the Quantum Hamilton equations" to solve the coupled FBSDE 
for a freely spinning particle with the gradient method or the
numerical approximation of the conditional expectation.
* Output: The results are written into .txt files
*/
#include <chrono>
#include <cmath>
#include "random"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <math.h>
#include <fstream>
#include <iterator>

using namespace std;

using T = double;

#include <stdio.h>
#include <stdlib.h>



enum method { euler, heun, runge };
enum problem {spin_f};
enum bsde_approach {u_conditional, q_conditional, q_gradient};
std::string problem_string[1] = {"spin_f"};
enum gauss_method {box_muller, standard, dunweg_paul};

void calculate_gradient_1d(T function[], int size, T grad[], T del_x) {
      // 4th order finite difference gradient
      grad[0] = -(25.0f*function[0] - 12.0f*(4.0f*function[1] - 3.0f*function[2]) - 16.0f*function[3] + 3.0f*function[4])/(12.0f*del_x);
      grad[1] = -(25.0f*function[1] - 12.0f*(4.0f*function[2] - 3.0f*function[3]) - 16.0f*function[4] + 3.0f*function[5])/(12.0f*del_x);
      grad[size-1] = (25.0f*function[size-1] - 12.0f*(4.0f*function[size-2] - 3.0f*function[size-3]) - 16.0f*function[size-4] + 3.0f*function[size-5])/(12.0f*del_x);
      grad[size-2] = (25.0f*function[size-2] - 12.0f*(4.0f*function[size-3] - 3.0f*function[size-4]) - 16.0f*function[size-5] + 3.0f*function[size-6])/(12.0f*del_x);
      for (size_t ind = 2; ind < size-2; ind++) {
              grad[ind] = (function[ind-2] - 8.0f*(function[ind-1] - function[ind+1]) - function[ind+2])/(12.0f * del_x);
      }
      return;
}

// get current bin oisition of position 1D
unsigned int get_x_bin_pos_1d (T x, unsigned int intervals, T dx, T x_low) {
  return (unsigned) int((x - x_low)/dx);
}

// get current bin position
unsigned int get_x_bin_pos (T x[], unsigned int intervals[], T dx[], T x_low[], uint8_t dimension) {
  unsigned int result = 0;
  for (size_t i = 0; i < dimension; i++) result += i*intervals[i]+ get_x_bin_pos_1d (x[i],intervals[i],dx[i],x_low[i]);
  return result;
}

// calculate q conditional for forward problem 1D
void calculate_q_conditional(T velt[], T x[], T dw[], T q_cond[], T q_cond_prev_iter_vec[], T q_cond_func[], unsigned int samples, unsigned int intervals, T dx, T x_low, T dt) {
  unsigned int x_bin = 0;
  unsigned int count[intervals+1];
  for (size_t ind = 0; ind < intervals; ind++) { q_cond[ind] = 0.0; count[ind] = 0;}
  for (size_t i = 0; i < samples; i++) {
        x_bin = get_x_bin_pos_1d(x[i+samples],intervals,dx,x_low);
        if ( x_bin < intervals && x_bin >= 0) {
              count[x_bin]++;
              q_cond[x_bin] += dw[i+samples] * velt[i];
        }
  }
  for (size_t ind = 0; ind < intervals; ind++) {
        if (count[ind]>0) q_cond[ind] /= (count[ind] * dt);
        else q_cond[ind] = q_cond_prev_iter_vec[ind];
        q_cond_func[ind]+=q_cond[ind];
  }
  return;
}
// check wheter theta angle is still in allowed segment
#define PI           3.14159265358979323846
bool check_theta(T val) {
    if (val > 0. && val < PI) {
        return true;
    }
    else {return false;}
}


int main (int argc, char* argv[])
{
                srand(time(NULL));
                std::cout<<"##############################Start##############################\n";

                // problem:
                problem p = problem::spin_f;
                // method of integration:
                method m = method::euler;
                // method for gaussian distributed random numbers:
                gauss_method g_m = gauss_method::box_muller;
                // approach for solveing the bsde
                bsde_approach bsde_a = bsde_approach::u_conditional;
                // T timestep:
                const T timestep = 0.001;
                // T mean of distribution:
                const T mean = 0.0;
                // int steps:
                const unsigned int steps = 1000;
                const unsigned int steps_p = steps+1; // steps + 1
                // int number of particles:
                const unsigned int particles = 100000;
                const T accuracy_bounds = 0.00001;
                // T lower bound of histogram:
                const T box_l = 0.0 + accuracy_bounds;
                // T upper bound of histogram:
                const T box_u = PI - accuracy_bounds;
                // int subInt_pos: number of intervals in spatial space between [box_l,box_u]
                const unsigned int subInt_pos = 200;
                // int iterations: number of iterations if FBSDEs are included
                const unsigned int iterations = 5;
                // print only 5 files
                unsigned int number_of_output_files = 5; number_of_output_files = iterations >= number_of_output_files ? int(iterations / number_of_output_files) : 1;
                // int sim_n: number of simulations
                const unsigned int sim_n = 1;

                unsigned int intervals = subInt_pos;

                // for random generators
                // define some constants
                const T dun_1 = (mean - sqrt (3.0 * timestep));
                const T dun_2 = (mean + sqrt (3.0 * timestep));
                const T sigma = sqrt( timestep );
                const T sq2 = sqrt( 2.0 );
                const T delta_x = (box_u - box_l)/(double(subInt_pos)); // width of the steps
                const T range = box_u - box_l;

                // double_well_fb constants
                // V(x) = V0/(x0^4) (x^2 - x0^2)^2
                // V0
                const T V0 = 2.0;//V00*0.5;
                // x0
                const T x0 = 2.0;//aa*0.5;
                // V'(x) = 4*V0/(x0^4) * x * (x^2 - x0^2) = const_dw1 * x (x^2 - const_dw2)
                const T const_dw1 = 4.0*V0/(pow(x0,4)), const_dw2 = x0*x0;

                // poeschl_teller constants
                const T lam = 1.0;
                const T lam0 = lam*(lam+1.0)*0.5;

                // spin  constants
                const T vchi0 = 0.0; // \Omega^v_chi
                const T vphi0 = 1.0; // \Omega^v_phi
                const T inertia = 1.0;
                const T inv_inertia = 1.0 / inertia;
                const T sq_inv_inertia = sqrt(inv_inertia);

                T r_mean;
                std::function<T (T)> startf;
                std::function<T (T)> vf;

                switch (p) {
                case spin_f:
                        startf=[range](T a){
                                    return 0.05+(PI-0.05)*a;
                            };
                        vf=[vphi0, vchi0, inv_inertia](T a){
                                    return -inv_inertia/(sin(a)*sin(a)*sin(a))*(-cos(a)* (vchi0 * vchi0 + vphi0 * vphi0) - (1+cos(a)*cos(a)) * vchi0 * vphi0);
                            };
                        break;
                default:
                        std::cout << "Invalid selection of method for gradient of the potential" << '\n';
                        break;
                }

                std::function<T (T)> random;
                random=[mean,sigma](int a){
                  // really slow !
                                std::mt19937 generator;
                                std::normal_distribution<T> distribution(mean,sigma);
                                return distribution(generator);
                        };

                switch (g_m)
                {
                case standard:
                        break;
                case  dunweg_paul:
                        random=[dun_1,dun_2](int a){
                                        T result = 0.0; T x_rand1;
                                        x_rand1 = 6.0*(static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)));
                                        if (x_rand1 < 1.0) {
                                                result = dun_1;
                                        } else if ( x_rand1 >= 5.0) {
                                                result = dun_2;
                                        }
                                        return result;
                                };
                        break;
                case box_muller:
                      random=[mean,sigma](int a){
                              static T n2 = 0.0;
                              static int n2_cached = 0;
                              if (!n2_cached)
                              {
                                  T x, y, r;
                                  do
                                  {
                                      x = 2.0*rand()/RAND_MAX - 1;
                                      y = 2.0*rand()/RAND_MAX - 1;
                                      r = x*x + y*y;
                                  }
                                  while (r == 0.0 || r > 1.0);
                                  {
                                      T d = sqrt(-2.0*log(r)/r);
                                      T n1 = x*d;
                                      n2 = y*d;
                                      T result = n1*sigma + mean;
                                      n2_cached = 1;
                                      return result;
                                  }
                              }
                              else
                              {
                                  n2_cached = 0;
                                  return n2*sigma + mean;
                              }
                      };
                      break;
                default:
                        //std::cout << "Invalid Selection of random number method" << '\n';
                        break;
                }

                        auto begin = std::chrono::system_clock::now();

                        T* path_iter= new T [particles*2];     // save all paths in array
                        T* path_iter_dW= new T [particles*2];          // save all dWs in array
                        T* stop_step= new T [particles];

                        // phi and chi angles
                        T* phi_= new T [particles*2];     // save all paths in array
                        T* dw_phi= new T [particles*2];          // save all dWs in array
                        T* stop_step_phi= new T [particles];
                        T* chi_= new T [particles*2];     // save all paths in array
                        T* dw_chi= new T [particles*2];          // save all dWs in array
                        T* stop_step_chi= new T [particles];


                        T* u_t_= new T [particles*2];         // osmotic velocity with respect to the time [particles][dimension]
                        T u_step_prev[intervals+1];         // osmotic_step_func from previous iteration plus under-/oveflow
                        T q_step_prev_iter_vec[intervals+1];
                        T u_step_func[intervals+1];         // osmotic_step_function[which Interval][dimension]
                        T q_step_func_vec[intervals+1];

                        T osmotic_step;     // osmotic velocity functional
                        T probability, probability_0;
                        T subInt_lb;
                        T conditional_step_func[intervals];

                        T osmotic_sim[intervals];     // osmotic velocity functional
                        T probability_sim[intervals], subInt_sim[intervals];

                        unsigned int cnt[intervals+1];

                        subInt_sim[0] = box_l + 0.5 * delta_x;
                        for (size_t in = 1; in < intervals; in++) {
                                subInt_sim[in] = subInt_sim[in-1]+delta_x;
                        }

                        int sum_osm[intervals+1];
                        T path_temp;
                        //T osmotic_temp_vec;
                        unsigned int int_pos_act, int_pos_act_next;

                        for (size_t sim_i = 0; sim_i < sim_n; sim_i++) {
                                        /************************************************************************/
                                        // Setting starting values for osmotic step function
                                        for (size_t ind = 0; ind < intervals+1; ind++) {
                                                sum_osm[ind] = 0.0;
                                                // does not work for hydrogen for example
                                                switch (p) {
                                                  default:
                                                        //u_step_prev[ind]=vf(box_l+delta_x*(ind+0.5));
                                                        // just for testing poeschl-teller
                                                        u_step_prev[ind]=0.0;
                                                        break;
                                                }
                                                q_step_prev_iter_vec[ind]=0.0;
                                                u_step_func[ind]=0.0;
                                                q_step_func_vec[ind]=0.0;
                                        } // END setting starting values for osmotic step function

                                        for (size_t iter = 0; iter < iterations; iter++) {
                                          int cou = 0;
                                                std::cout << "\r~>start iteration " << iter << "...                                                      " << std::flush;

                                                for (size_t j = 0; j < particles; j++) {
                                                        stop_step[j] = steps;
                                                        /******************** starting points ***************************/
                                                        path_iter[j]=startf(static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)));
                                                        //phi, chi
                                                        chi_[j]=2*PI*(static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)));
                                                        phi_[j]=2*PI*(static_cast <double> (rand()) /( static_cast <double> (RAND_MAX)));
                                                }
                                                /**************************** determine gradient of u prev iter******/
                                                T q_step[intervals];
                                                //calculate_gradient_1d(u_step_prev, intervals, q_step, delta_x);
                                                switch (bsde_a) {
                                                  case q_gradient:
                                                        calculate_gradient_1d(u_step_prev, intervals, q_step, delta_x);
                                                        break;
                                                  default:
                                                        break;
                                                }
                                                /********************************************************************/
/*************************************************** forward integration ***************/
                                                        for (size_t steps_j = 0; steps_j < steps-1; steps_j++) // START steps_j loop
                                                        /*******************************************************************************/
                                                        {
                                                          // gradient from previous iteration
                                                          if (!steps_j) {
                                                                  switch (bsde_a) {
                                                                    case q_gradient:
                                                                          calculate_gradient_1d(u_step_prev, intervals, q_step, delta_x);
                                                                          if (iter == iterations - 1) {
                                                                          for (size_t ind = 0; ind < intervals; ind++) std::cout << q_step[ind] << '\n';
                                                                          }
                                                                          break;
                                                                    default:
                                                                          break;
                                                                  }
                                                                  /********************************************************************/
                                                          }
                                                          // Calculate conditional_step[interval] for all intervals
                                                          T conditional_step[intervals];
                                                          switch (bsde_a) {
                                                                case q_gradient:
                                                                      for (size_t ind = 0; ind < intervals; ind++) {
                                                                              conditional_step[ind] = q_step[ind];
                                                                              conditional_step_func[ind]+=conditional_step[ind];
                                                                      }
                                                                      break;
                                                                case q_conditional:
                                                                      calculate_q_conditional(u_t_, path_iter, path_iter_dW, conditional_step, q_step_prev_iter_vec, q_step_func_vec, particles, intervals, delta_x, box_l, timestep);
                                                                      if (steps_j == steps-2 && iter == iterations -1) for (size_t ind = 0; ind < intervals; ind++) std::cout << ind << " " << conditional_step[ind] << '\n';
                                                                      break;
                                                                case u_conditional:
                                                                      for (size_t ind = 0; ind < intervals; ind++) {conditional_step[ind] = 0.0; cnt[ind] = 0;}
                                                                      // loop over all entries in interval_pos
                                                                      for (size_t jj = 0; jj < particles; jj++) {
                                                                            int_pos_act = get_x_bin_pos_1d(path_iter[jj+particles],intervals,delta_x,box_l);
                                                                            if ( int_pos_act < intervals && int_pos_act >= 0) {
                                                                                  conditional_step[int_pos_act] += u_t_[jj] + timestep * vf(path_iter[jj+particles]) ;
                                                                                  cnt[int_pos_act]++;// interval_pos[ind][steps_i][entry] = particle
                                                                            }
                                                                      }
                                                                      for (size_t ind = 0; ind < intervals; ind++) {
                                                                            if (cnt[ind]>0) conditional_step[ind] /= (cnt[ind]);
                                                                            else conditional_step[ind] = u_step_prev[ind];
                                                                            conditional_step_func[ind]+=conditional_step[ind];
                                                                      }
                                                                      if (steps_j == steps-2 && iter == iterations -1) for (size_t ind = 0; ind < intervals; ind++) std::cout << box_l+(0.5+ind)*delta_x << " " << conditional_step[ind] << '\n';
                                                                      break;
                                                                default:
                                                                      break;
                                                          }
                                                                for (size_t j = 0; j < particles; j++) { // START particle loop
                                                                      if (stop_step[j] > steps_j) {
                                                                            path_temp = path_iter[j];
                                                                            int_pos_act = (unsigned) int((path_temp - box_l)/delta_x );
                                                                            if (!steps_j) {
                                                                                    path_iter_dW[j] = random(0);
                                                                                    path_iter_dW[j+particles] = random(0);

                                                                                    // chi phi 
                                                                                             dw_chi[j] = random(0);
                                                                                             dw_phi[j] = random(0);   

                                                                                    // u^iter_new (start pos) = u^iter_old (start pos)
                                                                                    u_t_[j]  = u_step_prev[int_pos_act];
                                                                            }
                                                                            if ( !(int_pos_act < intervals) ) {
                                                                                    stop_step[j] = steps_j-1;                       // particle won't be included from step ii on in later calculation of osmotic velocity
                                                                                    break;
                                                                            }
                                                                            //if (j == 250 && (iter % number_of_output_files) == 1) std::cout << "x " << path_iter[j] << " step " << steps_j << '\n';
                                                                            // new position
                                                                            /**************************************************************/
                                                                            switch (p)
                                                                            {
                                                                            case spin_f:    
                                                                                path_iter[j+particles] = path_temp - inv_inertia * timestep * (u_t_[j] - 0.5*cos(path_temp)/sin(path_temp)) 
                                                                                                         + sq_inv_inertia * (cos(phi_[j]) * path_iter_dW[j] + sin(phi_[j]) * dw_phi[j]);
                                                                                phi_[j] += 0.5 * inv_inertia * timestep / (sin(path_iter[j]) * sin(path_iter[j])) * (1.0-cos(path_iter[j])) 
                                                                                           + sq_inv_inertia/sin(path_iter[j]) * (cos(path_iter[j])*cos(phi_[j])*dw_phi[j] - cos(path_iter[j])*sin(phi_[j])*path_iter_dW[j] + sin(path_iter[j])*dw_chi[j]);
                                                                                chi_[j] += 0.5 * inv_inertia * timestep / (sin(path_iter[j]) * sin(path_iter[j])) * (1.0-cos(path_iter[j])) 
                                                                                           + sq_inv_inertia * (sin(phi_[j])*path_iter_dW[j]-cos(phi_[j])*dw_phi[j]);

                                                                            default:
                                                                                    // going the u(x) way in integrating
                                                                                    //path_iter[j+particles] = path_temp + u_step_prev[int_pos_act] * timestep + path_iter_dW[j];
                                                                                    // use u(t) for integrating
                                                                                    path_iter[j+particles] = path_temp - u_t_[j] * timestep + path_iter_dW[j];
                                                                                    break;
                                                                            }
      /***************************************************START forward integration velocity***************/
                                                                            path_temp = path_iter[j];
                                                                            int_pos_act = get_x_bin_pos_1d(path_iter[j],intervals,delta_x,box_l);
                                                                            int_pos_act_next = get_x_bin_pos_1d(path_iter[j+particles],intervals,delta_x,box_l);
                                                                            switch (bsde_a) {
                                                                                case q_gradient:
                                                                                        u_t_[j+particles] = u_t_[j] + timestep * vf(path_temp) + sq_inv_inertia * q_step[int_pos_act] * random(0);//1.0/(cosh(path_iter[j])*cosh(path_iter[j])) * path_iter_dW[j];
                                                                                        break;
                                                                                case q_conditional:
                                                                                        u_t_[j+particles] = u_t_[j] + timestep * vf(path_temp) + sq_inv_inertia * conditional_step[int_pos_act] * path_iter_dW[j];
                                                                                        break;
                                                                                case u_conditional:
                                                                                        if (int_pos_act_next >= 0 && int_pos_act_next < intervals) u_t_[j+particles] = conditional_step[int_pos_act_next];// - timestep * vf(path_iter[j]);// + q_step[(unsigned) int((path_iter[j+particles] - box_l)/delta_x )];
                                                                                        break;
                                                                                default:
                                                                                        break;
                                                                            }
                                                                            switch (p) {
                                                                                case spin_f:
                                                                                        u_t_[j+particles] /= (1 - 0.5*inv_inertia * timestep/ (sin(path_temp)*sin(path_temp)));
                                                                                        break;
                                                                                default:
                                                                                        break;
                                                                            }
                                                                                // mean value of iterated osmotic_step_function:
                                                                            u_step_func[int_pos_act] += u_t_[j];
                                                                            q_step_func_vec[int_pos_act] += conditional_step[int_pos_act];
                                                                            sum_osm[int_pos_act] += 1;
                                                                            // übertrag
                                                                            u_t_[j] = u_t_[j+particles];
                                                                            path_iter_dW[j] = path_iter_dW[j+particles];
                                                                            //path_iter_dW[j] = random(0);
                                                                            path_iter_dW[j+particles] = random(0);
                                                                            path_iter[j] = path_iter[j+particles];
                                                                    } // stop_step check bracket
                                                                    else {
                                                                      cou++;
                                                                      if (steps_j == steps-2 && j == particles -2) { std::cout << "iter " << iter << " part " << j << " stop " << stop_step[j] << " c " << cou << '\n'; }
                                                                    }
                                                              } // END particle loop

                                                        } // END steps_j loop
                                                        std::ostringstream txt_name;
                                                        txt_name << "pltdata/"<<iter<<"_"<<problem_string[p]<<".txt";
                                                        std::string copyOfStr = txt_name.str();
                                                        const char* txt_name_ = copyOfStr.c_str();
                                                        ofstream myfile;
                                                        bool c_temp = (!(iter % number_of_output_files) or (iter == (iterations-1)));
                                                        if (c_temp) myfile.open (txt_name_);
                                                          /* index for last pos */
                                                          subInt_lb = box_l + 0.5 * delta_x;
                                                          probability_0 = 0;
                                                          T probability_norm = 0.0;
                                                          r_mean = 0.0;
                                                         // cartesian coordinates
                                                                  for (size_t ind = 0; ind < intervals; ind++) {
                                                                          if (!sum_osm[ind]) {
                                                                                  osmotic_step = u_step_prev[ind];
                                                                          }
                                                                          else {
                                                                                  u_step_prev[ind] = u_step_func[ind] / float(sum_osm[ind]);
                                                                                  q_step_prev_iter_vec[ind] = q_step_func_vec[ind] / float(sum_osm[ind]);
                                                                                  osmotic_step = u_step_prev[ind];
                                                                          }
                                                                          probability_0 += 2.0 * osmotic_step * delta_x; probability_norm += exp (probability_0) * delta_x;
                                                                  }
                                                                  probability_0 = 0;
                                                                  for (size_t ind = 0; ind < intervals; ind++) {
                                                                          if (!sum_osm[ind]) {
                                                                                  //      u_step_prev[ind] = 0;
                                                                                  osmotic_step = u_step_prev[ind];
                                                                          }
                                                                          else {
                                                                                  u_step_prev[ind] = u_step_func[ind] / float(sum_osm[ind]);
                                                                                  q_step_prev_iter_vec[ind] = q_step_func_vec[ind] / float(sum_osm[ind]);
                                                                                  osmotic_step = u_step_prev[ind];
                                                                          }
                                                                          probability_0 += 2.0 * osmotic_step * delta_x;
                                                                          std::cout << osmotic_step << " " << probability_0 << "\n";
                                                                          probability = exp (probability_0)/probability_norm;
                                                                                  osmotic_sim[ind] = osmotic_step;
                                                                                  probability_sim[ind] = probability;
                                                                          if (c_temp) myfile<<subInt_lb<<" "<<probability<<" "<<osmotic_step<<" "<<q_step_prev_iter_vec[ind]<<" "<<ind<<" "<<"\n";
                                                                          subInt_lb += delta_x;
                                                                  }
                                                          if (c_temp) myfile.close();
                                                          for (size_t ind = 0; ind < intervals; ind++) {sum_osm[ind] = 0; u_step_func[ind]=0.0; q_step_func_vec[ind] = 0.0;}
                                                /**************** END new particle paths *******************/
/***************************************************END forward integration ***************/


                                        } // END iteration loop
                        }// END FOR sim_n
                        std::cout << "\nr " << r_mean << '\n';
                        auto end = std::chrono::system_clock::now();
                        std::cout <<"\n~>Duration: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count()/1000.0f << " s\n" << std::endl;
                 //}} // V00 and aa loops
                 delete[] path_iter;     // save all paths in array
                 delete[] path_iter_dW;          // save all dWs in array
                 delete[] stop_step;
                 delete[] u_t_;
                std::cout << "###################################END#####################################" << '\n';

}
