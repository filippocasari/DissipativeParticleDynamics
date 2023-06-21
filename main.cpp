//
// Created by Filippo Casari on 28/04/23.
//
#include <iostream>
#include <thread>
#include <sstream>
#include <vector>
#include "Particle.h"
#include <unordered_set>
#include "matplotlib-cpp/matplotlibcpp.h"
#include <omp.h>
#include <map>
#include <random>
#define ITER_PLOT 10
using namespace std;
namespace plt = matplotlibcpp;


double compute_init_kinetic_energy(const std::vector<Particle> &particles) {
    double kinetic_energy = 0.0;
    for (auto &particle: particles) {
        kinetic_energy += particle.k_en;
    }
    return kinetic_energy;
}

struct PairHash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};




int main(int argc, char *argv[]) {
    //plt::detail::_interpreter::kill();
    int test =1; // 0 = no walls, just fluid, 1 = walls & chain molecules, 2 = fixed walls & ring molecules
    assert(test<3 && test>=0);
    double body_force =0.3;
    double L = 15;
    double rc = 1.0;
    int ro = 4;
    double sigma = 1.0;
    double KS =100;
    double RS;
    if(test==1){
        RS = 0.1;
    }else{
        RS = 0.3;
    }

    double gamma = 4.5;
    int N = ro * (int )L *(int) L;
    double m = 1.0;

    auto dt = static_cast<double>(std::strtod(argv[1], nullptr));
    std::vector<Particle> particle_list(N);


    if(test==2){
        particle_list.resize(N+10*9);

        double lower_bound = 1;
        double upper_bound = L-2;
        std::uniform_real_distribution<double> unif_chain(lower_bound,upper_bound);
        std::default_random_engine re_chain;

        int num_particles_in_molecules = 1000; //because 900 molecules of fluid and of walls
        for(int i=0; i<10; i++){
            double center_x = unif_chain(re_chain);
            double center_y = unif_chain(re_chain);
            assert(center_x>=0 && center_x<=L-1);
            assert(center_y>=0 && center_y<=L-1);
            //cout<< "center "<<center_x<<" "<<center_y<<endl;
            double delta_angle = 2*M_PI/9;
            for(int a=0; a<9;a++){
                double x = center_x + 1.0/ro * cos(a*delta_angle);
                double y = center_y + 1.0/ro * sin(a*delta_angle);
                int xcell = static_cast<int>( x / rc);
                int ycell = static_cast<int>( y / rc);
                //cout<<"chain "<<i<<" particle type= "<<"A, "<<xcell<<" "<<ycell<<endl;
                Particle p;
                if(a==0){
                    p= Particle(x,y,0,0,1,L,xcell,ycell,num_particles_in_molecules+a,"A", test, num_particles_in_molecules+8, num_particles_in_molecules+a+1, i);
                }else if(a==8){
                    p = Particle(x,y,0,0,1,L,xcell,ycell,num_particles_in_molecules+a,"A", test, num_particles_in_molecules+a-1, num_particles_in_molecules, i);
                }else{
                    p = Particle(x,y,0,0,1,L,xcell,ycell,num_particles_in_molecules+a,"A", test, num_particles_in_molecules+a-1, num_particles_in_molecules+a+1, i);
                }
                particle_list[N+a+i*9] = p;
            }
            num_particles_in_molecules+=9;



            }
    }


    if(test==1){
        particle_list.resize(N+42*7);

        double lower_bound = 1;
        double upper_bound = L-2;
        std::uniform_real_distribution<double> unif_chain(lower_bound,upper_bound);
        std::default_random_engine re_chain;

        int num_particles_in_molecules = 1000; //because 900 molecules of fluid and of walls
        for(int i=0; i<42; i++){
            double center_x = unif_chain(re_chain);
            double center_y = unif_chain(re_chain);
            assert(center_x>=0 && center_x<=L-1);
            assert(center_y>=0 && center_y<=L-1);
            cout<< "center "<<center_x<<" "<<center_y<<endl;
            for(int a=0; a<2;a++){
                int xcell = static_cast<int>( center_x+double(a)*1/ro / rc);
                int ycell = static_cast<int>( center_y / rc);
                //cout<<"chain "<<i<<" particle type= "<<"A, "<<xcell<<" "<<ycell<<endl;
                Particle p;
                if(a==0){
                    p= Particle(center_x+double(a)*1.0/ro,center_y,0,0,1,L,xcell,ycell,num_particles_in_molecules+a,"A", test, -1, num_particles_in_molecules+a+1, i);
                }else{
                    p = Particle(center_x+double(a)*1.0/ro,center_y,0,0,1,L,xcell,ycell,num_particles_in_molecules+a,"A", test, num_particles_in_molecules+a-1, num_particles_in_molecules+a+1, i);
                }


                particle_list[N+a+i*7] = p;
            }

            for(int b =2; b<7;b++){
                int xcell = static_cast<int>(center_x+double(b)*1/ro / rc);
                int ycell = static_cast<int>(center_y / rc);
                //cout<<"chain "<<i<<" particle type= "<<"B, "<<xcell<<" "<<ycell<<endl;
                Particle p;
                if(b ==6){
                    p = Particle(center_x+double(b)*1.0/ro,center_y,0,0,1,L,xcell,ycell,num_particles_in_molecules+b,"B", test, num_particles_in_molecules+b-1,-1,i);
                }else{
                    p = Particle(center_x+double(b)*1.0/ro,center_y,0,0,1,L,xcell,ycell,num_particles_in_molecules+b,"B", test, num_particles_in_molecules+b-1,num_particles_in_molecules+b+1,i);
                }

                particle_list[N+b+i*7] = p;
            }
            num_particles_in_molecules+=7;


        }

    }
    if(argv[1] == nullptr){
        cerr<<"Not enough arguments"<<endl;
        return -1;
    }

    printf("Assignment 4\n");
    cout << "Arguments to pass: dt" << endl;
    //signal(SIGINT, interrupt_handler);
    //int N = static_cast<int>(std::strtol(argv[1], nullptr, 10));

    //cout << "N: " << N << endl;



    std::vector<double> kinetic_energy_list;

    std::vector<double> tot_energy_list;
    std::vector<double> temperatures;
    std::vector<double> momentum_list;
    map <pair<string , string>, double> a_coeff;
    if(test ==0){
        a_coeff ={{{"F", "F"}, 25}};
    }
    if(test==1){
        a_coeff = {{{"A", "A"}, 50},{{"A", "B"}, 25}, {{"A", "F"}, 25}, {{"A", "W"}, 200}, {{"B", "A"}, 25}, {{"B", "B"}, 1}, {{"B", "F"}, 300}, {{"B", "W"}, 200},{ {"F", "A"}, 25}, { {"F", "B"}, 300}, { {"F", "F"}, 25}, { {"F", "W"}, 200}, { {"W", "A"}, 200}, { {"W", "B"}, 200} , { {"W", "F"}, 200}, { {"W", "W"}, 0}};
    }else{
        a_coeff = {{{"A", "A"}, 50}, {{"A", "F"}, 25}, {{"A", "W"}, 200}, {{"F", "A"}, 25}, {{"F", "F"}, 25}, {{"F", "W"}, 200}, {{"W", "A"}, 200}, {{"W", "F"}, 200}, {{"W", "W"}, 0} };
    }



    cout << "dt: " << dt << endl;


    int N_X = static_cast<int>(L / rc);
    int N_Y = static_cast<int>(L / rc);



    bool verbose = false;
    //double *velocities;

    //velocities = return_zero_velocities(N);
    //cout<< "velocities: " << velocities[0] << endl;


    printf("size vector particles: %lu\n", particle_list.size());


    cout<<"N = "<<N<<endl;
    std::uniform_real_distribution<double> unif(0.0,L);
    std::default_random_engine re;
    int N_fluid = 0;
    int N_walls = 0;
    for (int i = 0; i < N; i++) {
        //double alpha = 2.0 * 3.14 * rand();

        double x = unif(re);
        double y =unif(re);


        double vx = 0;
        double vy = 0;
        int xcell = static_cast<int>(x / rc);
        int ycell = static_cast<int>(y / rc);
        if((xcell ==0 || xcell==L-1) && test !=0){
            N_walls++;
        }else{
            N_fluid++;
        }
        if (verbose) {
            cout<< "particle: " << i << endl;
            cout << "xcell: " << xcell << " = " << x << " / " << rc << endl;
            cout << "ycell: " << ycell << " = " << y << " / " << rc << endl;
        }

        particle_list[i]=Particle(x, y, vx, vy, m, L, xcell, ycell, i , "F", test);

    }
    cout<< "particle list size: " << particle_list.size() << endl;

    std::unordered_set<std::pair<double, double>, PairHash> unique_positions;
    for (const auto &part: particle_list) {
        unique_positions.insert(std::make_pair(part.x, part.y));
    }
    if (unique_positions.size() == particle_list.size()) {
        std::cout << "All particles have unique coordinates" << std::endl;
    } else {
        std::cout << "Some particles have duplicate coordinates" << std::endl;
    }
    if (verbose) {
        for (int i = 0; i < (int) particle_list.size(); i++) {
            cout << "particle " << i << " position = " << particle_list[i].x << ", " << particle_list[i].y << endl;
            //cout << "particle " << i << " velocities = " << particle_list[i].vx << ", " << particle_list[i].vy << endl;
            cout << "particle " << i << " cell = " << particle_list[i].i << ", " << particle_list[i].j << endl;
        }
    }



    cout << "Initial kinetic energy: " << compute_init_kinetic_energy(particle_list) << endl;
    std::vector<std::vector<std::vector<Particle>>> cell_list(N_X, std::vector<std::vector<Particle>>(N_Y,
                                                                                                              std::vector<Particle>()));
    cout<<"size of cell list: "<<cell_list.size()<<endl;
    cout<<"size particle list "<<particle_list.size()<<endl;
    for (auto &i: particle_list) {
        //cout<<"particle "<<i.id<<" in cell "<<i.i<<", "<<i.j<<endl;
        cell_list[i.i][i.j].push_back(i);

    }
    cout<<"cell list created"<<endl;
    //assert(N_fluid+N_walls == (int) particle_list.size());

    vector<double> x_array_fluid = vector<double>(N_fluid);
    vector<double> y_array_fluid = vector<double>(N_fluid);
    vector<double> x_array_walls = vector<double>(N_walls);
    vector<double> y_array_walls = vector<double>(N_walls);
    vector<double> x_array_A ;
    vector<double> y_array_A;
    vector<double> x_array_B = vector<double>(42*5);
    vector<double> y_array_B = vector<double>(42*5);
    if(test==1){
        x_array_A= vector<double>(42*2);
        y_array_A = vector<double>(42*2);
    }else if(test==2){
        x_array_A= vector<double>(10*9);
        y_array_A = vector<double>(10*9);
    }




    for (int i = 0; i < (int) particle_list.size(); i++) {
        if(particle_list[i].type=="F"){
            x_array_fluid[i] = particle_list[i].x;
            y_array_fluid[i] = particle_list[i].y;
        }
        if(particle_list[i].type=="W"){
            x_array_walls[i] = particle_list[i].x;
            y_array_walls[i] = particle_list[i].y;
        }
        if(particle_list[i].type=="A"){
            x_array_A[i] = particle_list[i].x;
            y_array_A[i] = particle_list[i].y;
        }
        if(particle_list[i].type=="B"){
            x_array_B[i] = particle_list[i].x;
            y_array_B[i] = particle_list[i].y;
        }

    }
    cout<<"Num particles for Walls "<< N_walls<<endl;
    cout<<"Num particles for Fluid "<< N_fluid<<endl;
    cout<<"check Num particles for Walls "<< x_array_walls.size()<<endl;
    cout<<"check Num particles for Fluid "<< x_array_fluid.size()<<endl;
    cout<<"check Num particles for A "<< x_array_A.size()<<endl;
    if(test==1){
        cout<<"check Num particles for B "<< x_array_B.size()<<endl;
    }

    map<string, string> cmap_walls;
    cmap_walls.insert(make_pair("color", "gray"));
    map<string, string> cmap_A;
    cmap_A.insert(make_pair("color", "red"));
    map<string, string> cmap_B;
    cmap_B.insert(make_pair("color", "blue"));
    string color = cmap_walls["color"];
    cout<< " color of walls: "<<color<<endl;
    cout<< " color of A: "<<cmap_A["color"]<<endl;
    cout<< " color of B: "<<cmap_B["color"]<<endl;



    double T = 0.0;
    double T_dpd_thermostat = (sigma *sigma)/(2.0*gamma);
    cout << "T_dpd_thermostat: " << T_dpd_thermostat << endl;
    int iter = -1;

    auto start_time = std::chrono::high_resolution_clock::now();

    map <pair<int , int>, double> xi_map;

    normal_distribution<double> distribution(0.0,1.0);
    std::default_random_engine generator;
    for (int i=0; i< (int)particle_list.size()-1; i++) {
        for (int j=i+1; j< (int)particle_list.size(); j++) {
            assert(particle_list[i].id != particle_list[j].id);
            xi_map[make_pair(particle_list[i].id,particle_list[j].id)] = distribution(generator);


        }
    }
    cout<<"xi_map size: "<< xi_map.size()<<endl;
    //plt::clf();
    plt::figure();
    //plt::backend("TkAgg");
    plt::ion();

    plt::figure_size(800, 750);
    plt::title("Particle Methods with dt = "+to_string(dt)+" , T = "+to_string(T));



    while (iter < 20000) {
        //cout<<"starting the loop"<<endl;
        iter++;



        for (auto &part: particle_list) {
            int i_old = part.i;
            int j_old = part.j;
            double x_old = part.x;
            double y_old = part.y;

            part.update(dt, rc);

            int x_cell = part.i;
            int y_cell = part.j;

            cell_list[x_cell][y_cell].push_back(part);

            int in = -1;
            for (int i = 0; i < (int) cell_list[i_old][j_old].size(); i++) {
                if (cell_list[i_old][j_old][i].x == x_old && cell_list[i_old][j_old][i].y == y_old &&
                    cell_list[i_old][j_old][i].i == i_old && cell_list[i_old][j_old][i].j == j_old) {
                    in = i;
                    break;
                }
            }

            if (in != -1) {
                cell_list[i_old][j_old].erase(cell_list[i_old][j_old].begin() + in);

            } else {
                cout << "not found" << endl;
            }


        }



        for (auto &part: particle_list) {
            if(part.type =="W"){
                continue;
            }
            double x = part.x, y = part.y;
            int cell_x = part.i;
            int cell_y = part.j;
            double F_tot_x = 0.0;
            double F_tot_y = 0.0;
            if(test==2){
                F_tot_y+= body_force;
            }

            double F_x_C;
            double F_y_C ;
            double F_x_D = 0.0;
            double F_y_D = 0.0;
            double F_x_R ;
            double F_y_R ;
            double xi;
            double w_d;
            double d_v_x;
            double d_v_y;


            for (int i = -1; i < 2; i++) {

                for (int j = -1; j < 2; j++) {
                    int cell_x_new = (cell_x + i) % N_X;
                    int cell_y_new = (cell_y + j) % N_Y;
                    if (cell_x_new < 0) {
                        cell_x_new += N_X;
                    }

                    if (cell_y_new < 0) {
                        cell_y_new += N_Y;
                    }

                    //cout<< "num threads "<<omp_get_num_threads()<<endl;
                    for (auto &part2: cell_list[cell_x_new][cell_y_new]) {

                        double d_x = x - part2.x;
                        if (d_x > L / 2.0) {
                            d_x -= L;
                        } else if (d_x <= -L / 2.0) {
                            d_x += L;
                        }
                        double d_y = y - part2.y;
                        if (d_y > L / 2.0) {
                            d_y -= L;
                        } else if (d_y <= -L / 2.0) {
                            d_y += L;
                        }

                        double r = std::sqrt(d_x * d_x + d_y * d_y);


                        double w_r = 0.0;
                        double d_x_hat = d_x/r;
                        double d_y_hat = d_y/r;
                        if(r>0.0 ){
                            if (r <= rc) {
                                double a_i_j = a_coeff[{part.type, part2.type}];
                                F_x_C = a_i_j*(1.0 - r/rc)*d_x_hat;
                                F_y_C = a_i_j*(1.0 - r/rc)*d_y_hat;
                                w_r = 1.0 - r /rc;
                                w_d = w_r* w_r;
                                d_v_x = part.vx - part2.vx;
                                d_v_y = part.vy - part2.vy;

                                F_x_D = -gamma * w_d * (d_x_hat)* (d_x_hat*d_v_x+d_y_hat*d_v_y);
                                F_y_D = -gamma * w_d * (d_y_hat)* (d_x_hat*d_v_x+d_y_hat*d_v_y);

                                if (xi_map.find({part.id, part2.id}) == xi_map.end()) {
                                    // not found
                                    xi = xi_map[{part2.id, part.id}];

                                } else {
                                    // found
                                    xi = xi_map[{part.id, part2.id}];
                                }

                                F_x_R = sigma * w_r * xi * d_x_hat;
                                F_y_R = sigma * w_r * xi * d_y_hat;


                                F_tot_x+= (F_x_C) + (F_x_D) + (F_x_R/(sqrt(dt)));
                                F_tot_y+= F_y_C + F_y_D + F_y_R/(sqrt(dt));

                            }
                            if(((part2.type == "A") || part2.type == "B") && part.bond_id == part2.bond_id && (part.left_p == part2.id || part.right_p == part2.id)){
                                double F_x_S = KS *(1.0-r/RS)*d_x_hat;
                                double F_y_S = KS *(1.0-r/RS)*d_y_hat;
                                F_tot_x+= F_x_S;
                                F_tot_y+= F_y_S;
                            }

                        }


                    }
                }

            }
            part.update_vel_acc(dt, F_tot_x, F_tot_y);

        }

        double kinetic_energy = 0;
        double momentum;
        double sum_vx = 0;
        double sum_vy = 0;
        //cout<<"before counting particles within the loop"<<endl;
        int counter_wall =0;
        int counter_fluid =0;
        int counter_A =0;
        int counter_B =0;

        for (int i = 0; i < (int)particle_list.size(); i++) {
            //cout<<"particle type "<<particle_list[i].type<<endl;
            if(particle_list[i].type =="F"){

                x_array_fluid[counter_fluid] = particle_list[i].x;
                y_array_fluid[counter_fluid] = particle_list[i].y;
                counter_fluid++;
            }else if(particle_list[i].type =="W"){
                x_array_walls[counter_wall] = particle_list[i].x;
                y_array_walls[counter_wall] = particle_list[i].y;
                counter_wall++;
            }else if (particle_list[i].type =="A"){
                x_array_A[counter_A] = particle_list[i].x;
                y_array_A[counter_A] = particle_list[i].y;
                counter_A++;
            }else if (particle_list[i].type =="B"){
                x_array_B[counter_B] = particle_list[i].x;
                y_array_B[counter_B] = particle_list[i].y;
                counter_B++;
            }


            kinetic_energy += particle_list[i].k_en;

            sum_vx += particle_list[i].vx;
            sum_vy += particle_list[i].vy;

        }

        kinetic_energy_list.push_back(kinetic_energy);
        T = 2.0 * kinetic_energy / (3.0 * (double)N);
        temperatures.push_back(T);

        momentum = sqrt(sum_vx * sum_vx + sum_vy * sum_vy) * m;
        momentum_list.push_back(momentum);
        //verbose =true;
        if(verbose){cout<<"iter: "<< iter<<endl;}

        if (iter % ITER_PLOT == 0) {


            //plt::clf();
            if(test==0){
                plt::subplot(3, 1, 1);
                plt::cla();
                plt::xlim(0.0, (double) L);
                plt::ylim(0.0, (double) L);
                plt::grid(true);
                plt::scatter(x_array_fluid, y_array_fluid, 20);

                plt::subplot(3, 1, 2);

                plt::named_plot("temperature", temperatures, "r");
                if (iter == 0) { plt::legend({{"loc", "upper right"}}); }
                plt::grid(true);
                plt::xlabel("Iterations");
                plt::subplot(3, 1, 3);
                plt::named_plot("momentum", momentum_list, "blue");
                if (iter == 0) { plt::legend({{"loc", "upper right"}}); }
                plt::grid(true);
            }else{
                plt::subplot(2, 1, 1);
                plt::cla();
                plt::xlim(0.0, (double) L);
                plt::ylim(0.0, (double) L);
                plt::grid(true);
                plt::scatter(x_array_fluid, y_array_fluid, 20);

                plt::scatter(x_array_walls, y_array_walls, 20,cmap_walls );
                plt::scatter(x_array_A, y_array_A, 20,cmap_A );
                if( test ==1){
                    plt::scatter(x_array_B, y_array_B, 20,cmap_B );
                }
                plt::subplot(2, 1, 2);

                plt::named_plot("temperature", temperatures, "r");
                if (iter == 0) { plt::legend({{"loc", "upper right"}}); }
                plt::grid(true);
                plt::xlabel("Iterations");
            }

            //plt::subplot(2, 2, 2);
            //plt::named_plot("kinetic", kinetic_energy_list, "red");

            //plt::xlabel("Iterations");
            //
            //plt::grid(true);

            //plt::subplot(2, 2, 4);

            //if (iter == 0) { plt::legend({{"loc", "upper right"}}); }
            //plt::grid(true);
            //plt::xlabel("Iterations");
            //if (iter == 0) { plt::legend({{"loc", "upper right"}}); }
            plt::draw();

            plt::pause(0.001);

        }


    }
    auto end = std::chrono::high_resolution_clock::now(); // get end time
    auto duration = (end - start_time); // calculate duration

    cout <<endl<<"Execution time: " << (double)duration.count() / 1000000.0 << " seconds" << endl;

    plt::detail::_interpreter::kill();

    cout << "AS4 completed" << endl;
    return 0;
}



