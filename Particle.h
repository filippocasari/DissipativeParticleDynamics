//
// Created by Filippo Casari on 26.04.23.
//

#ifndef AS3PM_PARTICLE_H
#define AS3PM_PARTICLE_H


#include <cmath>
#include <iostream>
#include <utility>
class Particle
{
    // type 0 = fluid, 1=wall
public:
    double x;
    std::string type = "F";
    double y;
    double vx;
    double vy;
    double m;
    double L;
    double a_x = 0;
    double a_y = 0;
    double F_x=0 ;
    double F_y=0 ;
    int i;
    int j;
    int id;
    double vel_walls_y =5.0;
    int left_p = -1;
    int right_p = -1;
    int bond_id ;
    double k_en=0;
    Particle() = default;
    Particle(double x, double y, double vx, double vy, double m, double L, int i, int j, int id, std::string type, int test, int left_p=-1, int right_p=-1, int bond_id=-1)
    {
        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
        this->m = m;
        this->L = L;
        this->i = i;
        this->j = j;
        this->type = type;
        this->id = id;
        if(this->type=="F" && test!=0){
            if(((i ==0) || (i ==L-1))){
                this->type="W";
                if(this->i==0){
                    this->vel_walls_y = -this->vel_walls_y;

                }
                if(test == 2){
                    this->vel_walls_y = 0;
                }

                this->k_en = 0.5 * this->m * (this->vel_walls_y * this->vel_walls_y);
            }
            else {
                this->compute_k_energy();
            }

        }
        else{
            this->compute_k_energy();
        }
        if(this->type=="A" || this->type=="B"){
            if(left_p!=-1){
                this->left_p = left_p;
            }
            if(right_p!=-1){
                this->right_p = right_p;
            }
            this->bond_id = bond_id;
        }
        else{
            this->bond_id = -1;
        }


    }

    void update(double dt, double rc)
    {
        if(this->type == "W"){
            this->y += (this->vel_walls_y * dt);

        }else{
            this->y += (this->vy * dt) + this->a_y * dt * dt * 0.5;
            this->x += (this->vx * dt) + this->a_x * dt * dt * 0.5;
            if(this->x < 0){
                this->x += this->L;
            }else if(this->x > this->L){
                this->x -= this->L;
            }
        }


        if(this->y < 0){
            this->y += this->L;
        }else if(this->y > this->L){
            this->y -= this->L;
        }

        this->i = int(this->x / rc);
        this->j = int(this->y / rc);
        if(this->i<0){
            this->i += int(this->L/rc);
        }
        if(this->j<0){
            this->j += int(this->L/rc);
        }
    }
    void update_vel_acc(double dt,  double Fx, double Fy)
    {
        if(this->type=="F" || this->type=="A" || this->type=="B"){
            this->F_x = Fx;
            this->F_y = Fy;
            //double tmp_vel = sqrt(this->vx*this->vx + this->vy*this->vy);
            double new_acc_x = this->F_x / this->m;
            double new_acc_y = this->F_y / this->m;
            this->vx += (this->a_x + new_acc_x) * dt * 0.5;
            this->vy += (this->a_y + new_acc_y) * dt * 0.5;
            this->a_x = new_acc_x;
            this->a_y = new_acc_y;
            this->compute_k_energy();
        }
        else if(type=="W"){
            this->k_en = 0.5 * this->m * (this->vel_walls_y * this->vel_walls_y);
        }

        //double new_vel = sqrt(this->vx*this->vx + this->vy*this->vy);
        //std::cout<<"difference velocity : "<< new_vel-tmp_vel<<std::endl;


    }
    inline void compute_k_energy(){
        this->k_en = 0.5 * this->m * (this->vx*this->vx + this->vy * this->vy);
    }

    bool operator==(const Particle& other) const {
        return this->x == other.x && this->y == other.y && this->vx == other.vx && this->vy == other.vy && this->m == other.m && this->L == other.L && this->i == other.i && this->j == other.j;
    }
    bool operator!=(const Particle& other) const {
        return !(*this == other);
    }


};


#endif
