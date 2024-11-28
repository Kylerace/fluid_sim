
#include <cmath>
#include <utility>
#include <stdio.h>
#include <windows.h>
#include <iostream>
#include <ctime>
#include <string>
#include <iterator>
#include <random>
#include <chrono>
#include <thread>
#include <fstream>
#include <sstream>
#include <filesystem>
//#include <minwindef.h>

//#include "window.h"
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
#include "dependencies/glad/glad.h"
#include "dependencies/GLFW/include/GLFW/glfw3.h"

#include <algorithm>
#include <optional>
#undef max

using namespace std;

typedef void (*GL_GENBUFFERS) (GLsizei, GLuint*);
//typedef void (*grav_pointer) (State*, fluid_type, fluid_type, int, int);
typedef float fluid_type;


bool sample_log(int expected_passes) {
    bool ret = false;
    int straw = (int)(rand() % expected_passes);
    //printf("\n\tstraw:%i", straw);
    if(straw < log10((double)expected_passes)) {
        ret = true;
    }
    return ret;
}

#define U_FIELD 0
#define V_FIELD 1
#define S_FIELD 2
#define T_FIELD 3


#define RENDER_MODE_PRESSURE 0
#define RENDER_MODE_VELOCITY 1
#define RENDER_MODE_TEMPERATURE 2

#define SCENE_TYPE_TANK 0
#define SCENE_TYPE_WIND_TUNNEL 1
#define SCENE_TYPE_PAINT 2

#define S_FIELD_FLUID 1.0
#define S_FIELD_SOLID 0.0

#define MIN_TIME_TO_LOG 3.0

#define DONT_AUTO_PAUSE -1.0
#define PAUSE_AFTER 10.0

#define WHITE_OUT_UNTOUCHED_CELLS FALSE

#define DEFAULT_DT 1.0 / 60.0
#define DEFAULT_VISCOSITY 0.0
#define DEFAULT_RESOLUTION 1000
#define DEFAULT_SCENARIO SCENE_TYPE_WIND_TUNNEL
#define DEFAULT_DENSITY 1000.0
#define DEFAULT_OVER_RELAXATION 1.6
#define DEFAULT_WIND_TUNNEL_IN_VEL 2.0
#define DEFAULT_DIVERGENCE_ITERATIONS 50
#define DEFAULT_DIVERGENCE_ITERATIONS_WIND_TUNNEL 100
#define DEFAULT_COARSE_GRID_BUDGET_DIVISOR 2
#define DEFAULT_COARSE_GRID_SIDE_LENGTH_POWER 3
#define DEFAULT_MINIMUM_POW_2_CELL_DIVISIBILITY 4

#define STARTING_DT DEFAULT_DT 
#define STARTING_VISCOSITY DEFAULT_VISCOSITY
#define STARTING_RESOLUTION 1000
#define STARTING_SCENARIO DEFAULT_SCENARIO
#define STARTING_DENSITY DEFAULT_DENSITY
#define STARTING_OVER_RELAXATION 1.8
#define STARTING_WIND_TUNNEL_IN_VEL DEFAULT_WIND_TUNNEL_IN_VEL
#define STARTING_DIVERGENCE_ITERATIONS DEFAULT_DIVERGENCE_ITERATIONS
#define STARTING_DIVERGENCE_ITERATIONS_WIND_TUNNEL 3000
#define STARTING_COARSE_GRID_BUDGET_DIVISOR 3
#define STARTING_COARSE_GRID_SIDE_LENGTH_POWER 1
#define STARTING_MINIMUM_POW_2_CELL_DIVISIBILITY 5

#define INITIAL_WINDOW_WIDTH 1920
#define INITIAL_WINDOW_HEIGHT 1080

bool print_perf = true;

int render_mode = RENDER_MODE_PRESSURE;

bool request_show_pressure = false;
bool request_show_smoke = false;
bool request_show_velocity = false;
bool request_show_streamlines = false;
bool request_show_temperature = false;
bool show_temperature_setting = false;
bool show_pressure_setting = false;

/// set by the glfw callback
bool change_mode_active = false;
string change_mode_buffer;

string viscosity_command = "VISCOSITY";
string divergence_iters_command = "DIVERGENCE";
string pulse_command = "PULSE"; //creates a pulse of pressure along the left side of the world with specified x vel
string grav_power_command = "GRAV";

#define COMMAND_VISCOSITY "VISCOSITY"
#define COMMAND_DIVERGENCE "DIVERGENCE"
#define COMMAND_PULSE "PULSE"
#define COMMAND_GRAV "GRAV"
#define COMMAND_PRINT_PERF "PRINT PERF"

vector<string> commands{COMMAND_VISCOSITY, COMMAND_DIVERGENCE, COMMAND_PULSE, COMMAND_GRAV, COMMAND_PRINT_PERF};

bool request_clear_pressure = false;

int request_delta_pressure_iterations = 0;

int request_scene_type = SCENE_TYPE_TANK;

bool paused = false;

struct Object {
    fluid_type x;
    fluid_type y;
    fluid_type radius;
};

struct State {
    fluid_type density;
    int* cells;
    int dimensions;
    int num_cells;
    /// cell side length
    fluid_type h;

    fluid_type* u_current; //u[our index] = velocity at our leftmost edge
    fluid_type* v_current; //v[our index] = velocity at our down edge
    fluid_type* u_next; //whats actually operated on in sim_step(), switch with current at the end
    fluid_type* v_next;

    fluid_type* pressure; //helps with enforcing 0 divergence
    fluid_type* solids; //the walls and the circle

    fluid_type* u_source;
    fluid_type* v_source;

    fluid_type* m_current; //smoke
    fluid_type* m_next;

    fluid_type* temperature_current; //independent from velocity
    fluid_type* temperature_next;

    int coarse_cells_x;
    int coarse_cells_y;
    int coarse_cells;
    int fine_cells_per_coarse;
    int coarse_grid_side_length_mult;
    fluid_type coarse_h;

    fluid_type* u_coarse;
    fluid_type* v_coarse;
    fluid_type* u_residuals;
    fluid_type* v_residuals;
    fluid_type* solids_coarse;

    /// reset to 0 every sim step, set to not 0 to force that pixel to a color
    int* debug_touch_buffer;
    /// idk what this is for
    fluid_type* debug_data_buffer;

    fluid_type viscosity;
    fluid_type over_relaxation;
    int num_divergence_iterations;
    fluid_type g_strength;
    void (*g_pointer) (State* s, fluid_type strength, fluid_type dt, int x, int y);
    int objects_len;
    int objects_capacity;
    Object* objects;
    bool clear_pressure;
};

void set_obstacle(State* s, fluid_type dt, int index, fluid_type x, fluid_type y, bool reset = false) {
    if(index >= s->objects_len) {
        return;
    }
    Object* object = &s->objects[index];

    fluid_type vx = 0.0;
    fluid_type vy = 0.0;

    if(!reset) {
        vx = (x - object->x) / dt;
        vy = (y - object->y) / dt;
    }

    object->x = x;
    object->y = y;

    int n = s->cells[0];
    fluid_type cd = sqrt(2) * s->h;
    fluid_type r = object->radius;

    for(int j = 1; j < s->cells[1] - 2; j++) {
        for(int i = 1; i < s->cells[0] - 2; i++) {
            s->solids[i + n * j] = S_FIELD_FLUID;

            fluid_type d = sqrt(pow(x - i, 2.0) + pow(y - j, 2.0));

            //if(dx * dx + dy * dy < r * r) {
            if(d < r) {
                //printf("\n%g, %g outside of object with radius %g lhs: %g rhs: %g", dx, dy, r, dx * dx + dy + dy, r * r);
                //int x = 1 / 0;
                s->solids[j * n + i] = S_FIELD_SOLID;

                s->m_current[j * n + i] = 1.0;

                s->u_current[j * n + i] = vx;
                s->u_current[j * n + i + 1] = vx;
                s->v_current[j * n + i] = vy;
                s->v_current[(j + 1) * n + i] = vy;
            }
        }
    }
}

void add_object(State* s, fluid_type x, fluid_type y, fluid_type radius) {
    if(s->objects_capacity == s->objects_len) {

        s->objects_capacity *= 2;
        Object* new_objects = new Object[s->objects_capacity];
        copy(s->objects, s->objects + s->objects_len - 1, new_objects);
        delete [] s->objects;
        s->objects = new_objects;
    }

    Object* object = &s->objects[s->objects_len];
    object->radius = radius;
    object->x = x;
    object->y = y;
    s->objects_len++;
    set_obstacle(s, 0.0, s->objects_len - 1, x, y, true);
}

template<typename T>
void zero_out(T* to_zero, int len) {
    for(int i = 0; i < len; i++) {
        to_zero[i] = (T)0.0;
    }
}

void vary_averaged_quantities(State* s, int x, int y) {

    double lower_bound = -1.0;
    double upper_bound = 1.0;

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();

    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine re(seed);

    double periods_per_meter[2] = {20.0, 10.0};
    double mean_temp = 350.0;
    const double pi = 3.14159265;

    double sin_x = 5.0 * sin(2 * pi * periods_per_meter[0] * (x * s->h)); 
    double sin_y = 3.0 * sin(2 * pi * periods_per_meter[1] * (y * s->h));

    s->temperature_current[x + y * s->cells[0]] = mean_temp + sin_x + sin_y + unif(re);
    s->temperature_next[x + y * s->cells[0]] = s->temperature_current[x + y * s->cells[0]];
}

State* new_state(fluid_type density, fluid_type viscosity, int num_x, int num_y, fluid_type h, fluid_type over_relaxation, fluid_type g_strength, void (*g_pointer) (State*, fluid_type, fluid_type, int, int)) {
    State* s = new State();
    s->density = density;
    s->viscosity = viscosity;

    num_x += 2;
    num_y += 2;

    //each dimension must be divisible by this
    int modulo = (int)pow(2, STARTING_MINIMUM_POW_2_CELL_DIVISIBILITY);

    while(num_x % modulo != 0) {
        num_x++;
    }
    while(num_y % modulo != 0) {
        num_y++;
    }

    s->cells = new int[2] {num_x, num_y};
    s->dimensions = 2;
    s->h = h;
    s->over_relaxation = over_relaxation;
    s->num_cells = s->cells[0] * s->cells[1];
    printf("\nnum_cells: %i, x: %i, y: %i", s->num_cells, s->cells[0], s->cells[1]);
    s->u_current = new fluid_type[s->num_cells];
    s->u_source = new fluid_type[s->num_cells];
    s->v_current = new fluid_type[s->num_cells];
    s->v_source = new fluid_type[s->num_cells];
    s->u_next = new fluid_type[s->num_cells];
    s->v_next = new fluid_type[s->num_cells];
    s->pressure = new fluid_type[s->num_cells];
    s->solids = new fluid_type[s->num_cells];
    s->m_current = new fluid_type[s->num_cells];
    s->m_next = new fluid_type[s->num_cells];

    s->temperature_current = new fluid_type[s->num_cells];
    s->temperature_next = new fluid_type[s->num_cells];

    int coarse_grid_side_length_power = STARTING_COARSE_GRID_SIDE_LENGTH_POWER;
    int coarse_grid_side_length_mult = (int)pow(2, coarse_grid_side_length_power);
    s->coarse_grid_side_length_mult = coarse_grid_side_length_mult;

    s->coarse_cells_x = num_x / coarse_grid_side_length_mult;
    s->coarse_cells_y = num_y / coarse_grid_side_length_mult;
    s->coarse_cells = s->coarse_cells_x * s->coarse_cells_y;
    s->fine_cells_per_coarse = (int)pow(coarse_grid_side_length_mult, 2);
    s->coarse_h = h * coarse_grid_side_length_mult;

    s->u_coarse = new fluid_type[s->coarse_cells];
    s->v_coarse = new fluid_type[s->coarse_cells];
    s->u_residuals = new fluid_type[s->coarse_cells];
    s->v_residuals = new fluid_type[s->coarse_cells];
    s->solids_coarse = new fluid_type[s->coarse_cells];

    s->debug_touch_buffer = new int[s->num_cells];
    s->debug_data_buffer = new fluid_type[s->num_cells];

    s->num_divergence_iterations = 100;
    s->g_strength = g_strength;
    s->g_pointer = g_pointer;
    s->objects = new Object[5];
    s->objects_capacity = 5;
    s->objects_len = 0;
    s->clear_pressure = false;

    int num = s->num_cells;
    zero_out(s->u_current, num);
    zero_out(s->u_source, num);
    zero_out(s->v_current, num);
    zero_out(s->v_source, num);
    zero_out(s->u_next, num);
    zero_out(s->v_next, num);
    zero_out(s->pressure, num);
    zero_out(s->solids, num);
    zero_out(s->m_next, num);
    zero_out(s->debug_touch_buffer, num);
    zero_out(s->debug_data_buffer, num);
    zero_out(s->u_coarse, s->coarse_cells);
    zero_out(s->v_coarse, s->coarse_cells);
    zero_out(s->u_residuals, s->coarse_cells);
    zero_out(s->v_residuals, s->coarse_cells);
    zero_out(s->solids_coarse, s->coarse_cells);


    for(int y = 0; y < s->cells[1]; y++) {
        for(int x = 0; x < s->cells[0]; x++) {
            vary_averaged_quantities(s, x, y);
        }
    }

    for(int i = 0; i < num; i++){
        s->m_current[i] = 1.0;
    }
    return s;
}


void g_standard(State* s, fluid_type strength, fluid_type dt, int x, int y) {
    //*v += -strength * dt;
    s->v_current[y * s->cells[0] + x] += strength * dt;
}


void g_right(State* s, fluid_type strength, fluid_type dt, int x, int y) {
    if(x < s->cells[0] / 2) {
        s->v_current[y * s->cells[0] + x] += (strength / 1.0) * dt;
    } else {
        s->v_current[x + s->cells[0] * y] += -(strength / 1.0) * dt;
    }
    //*v += -strength * dt;
    //s->u_current[x + s->cells[0] * y] += -strength * dt;
}

void g_center(State* s, fluid_type strength, fluid_type dt, int x, int y) {
    fluid_type theta = atan2(y - s->cells[1] / 2, x - s->cells[0] / 2);
    fluid_type magnitude = -strength / pow(max(sqrt(pow(x - s->cells[0] / 2, 2) + pow(y - s->cells[1] / 2, 2)), 0.05), 2);

    fluid_type x_mag = magnitude * cos(theta);
    fluid_type y_mag = magnitude * sin(theta);

    //printf("\n\t\t\tx, y: %i, %i x force: %g y force: %g", x, y, x_mag, y_mag);

    //fluid_type magnitude = -strength / max((int)pow(x - s->cells[0] / 2, 2) + (int)pow(y - s->cells[1] / 2,2), 1);
    s->u_current[x + s->cells[0] * y] += x_mag * dt;
    s->v_current[x + s->cells[0] * y] += y_mag * dt;
}

struct Scene {
    State* s;
    int scene_type;
    bool show_smoke;
    bool show_pressure;
    bool show_streamlines;
    bool show_velocity;
    bool show_temperature;

    public:

        Scene() {
            //State* s = new_state(1000.0, 1.0, 500, 500, 0.05, 1.9);
            int res = STARTING_RESOLUTION;
            s = new_state(STARTING_DENSITY, STARTING_VISCOSITY, res, res, 1.0 / res, STARTING_OVER_RELAXATION, 9.81, g_standard);
            scene_type = STARTING_SCENARIO;
            show_smoke = false;
            show_pressure = true;
            show_streamlines = false;
            show_velocity = false;
        }

        void set_request_globals() {
            request_scene_type = scene_type;
            request_show_pressure = show_pressure;
            request_show_smoke = show_smoke;
            request_show_streamlines = show_streamlines;
            request_show_velocity = show_velocity;
            request_show_temperature = show_temperature;
            show_pressure_setting = show_pressure;
            show_temperature_setting = show_temperature;
            request_clear_pressure = s->clear_pressure;
        }

        void setup_scene(int scene_type) {
            this->scene_type = scene_type;

            s->objects_len = 0;

            s->clear_pressure = true;

            s->num_divergence_iterations = STARTING_DIVERGENCE_ITERATIONS;
            for(int y = 0; y < s->cells[1]; y++) {
                for(int x = 0; x < s->cells[0]; x++) {
                    int index = x + s->cells[0] * y;
                    s->solids[index] = S_FIELD_FLUID;
                    s->pressure[index] = 0.0;
                    s->m_current[index] = 1.0;
                    s->m_next[index] = 0.0;
                    s->u_current[index] = 0.0;
                    s->u_next[index] = 0.0;
                    s->v_current[index] = 0.0;
                    s->v_next[index] = 0.0;
                    vary_averaged_quantities(s, x, y);
                }
            }
            
            if(scene_type == SCENE_TYPE_TANK) {
                for(int y = 0; y < s->cells[1]; y++) {
                    for(int x = 0; x < s->cells[0]; x++) {
                        fluid_type sol = S_FIELD_FLUID;
                        if(x == 0 || x == s->cells[0] - 1 || y == 0) {
                            sol = S_FIELD_SOLID;
                        }
                        s->solids[x + s->cells[0] * y] = sol;
                    }
                }
                s->num_divergence_iterations = 50;
                s->g_strength = -9.81;
                show_pressure = true;
                show_pressure_setting = true;
                show_smoke = false;
                show_streamlines = false;
                show_velocity = false;
                show_temperature = false;
                show_temperature_setting = false;
            }
            else if(scene_type == SCENE_TYPE_WIND_TUNNEL) {
                fluid_type in_velocity = STARTING_WIND_TUNNEL_IN_VEL;
                for(int y = 0; y < s->cells[1]; y++) {
                    for(int x = 0; x < s->cells[0]; x++) {
                        fluid_type sol = S_FIELD_FLUID;
                        if(x == 0 || y == s->cells[1] - 1 || y == 0) {
                            sol = S_FIELD_SOLID;
                        }
                        s->solids[x + s->cells[0] * y] = sol;

                        if(x == 1) {
                            s->u_current[x + s->cells[0] * y] = in_velocity;
                            //s->v_current[x + s->cells[0] * y] = in_velocity / 0.50;
                            //s->u_source[x + s->cells[0] * y] = in_velocity;
                        }
                        //if(x == (s->cells[0] - 3)) {
                        //    s->u_current[x + s->cells[0] * y] = 1.0 * in_velocity;
                        //}
                        //if(x == s->cells[0] - 1) {
                        //    s->u_current[x + s->cells[0] * y] = -in_velocity;
                        //}
                    }
                }

                s->num_divergence_iterations = STARTING_DIVERGENCE_ITERATIONS_WIND_TUNNEL;

                fluid_type pipe_h = 0.1 * s->cells[1];
                int min_j = floor(0.5 * s->cells[1] - 0.5 * pipe_h);
                int max_j = floor(0.5 * s->cells[1] + 0.5 * pipe_h);


                for(int y = min_j; y < max_j; y++) {
                    //for(int x = 0; x < s->cells[0]; x++) {
                    s->m_current[s->cells[0] * y] = 0.0;
                    s->m_current[s->cells[0] * y + s->cells[0] - 2] = 0.0;
                    //}
                }

                s->g_strength = 0.0;
                show_pressure = true;
                show_smoke = true;
                show_streamlines = false;
                show_velocity = false;
            } 
            else if(scene_type == SCENE_TYPE_PAINT) {
                s->g_strength = 0.0;
                s->over_relaxation = 1.0;
                show_pressure = false;
                show_smoke = true;
                show_streamlines = false;
                show_velocity = false;
            }

            printf("\n x at %i y at %i", (int)(s->cells[0] / 2.0), (int)(s->cells[1] / 2.0));
            add_object(s, (int)(s->cells[0] * 0.4), (int)(s->cells[1] / 2.0), s->cells[1] * 0.15);
            set_request_globals();

            zero_out(s->solids_coarse, s->coarse_cells);
            int n = s->cells[0];
            int coarse_cells_x = s->coarse_cells_x;
            int coarse_cells_y = s->coarse_cells_y;
            int coarse_side_multiplier = s->coarse_grid_side_length_mult;
            fluid_type* coarse_solids = s->solids_coarse;
            fluid_type* solids_buffer = s->solids;
            int cn = s->coarse_cells_y;
            for(int cy = 0; cy < coarse_cells_y; cy++) {
                for(int cx = 0; cx < coarse_cells_x; cx++) {
                    fluid_type u_residual = 0.0; //the sum of the divergence of each of the cells we're downsampled from
                    fluid_type v_residual = 0.0;
                    int cus = cx + coarse_cells_y * cy;
                    int cright = cus + 1;
                    int cleft = cus - 1;
                    int cup = cus + coarse_cells_y;
                    int cdown = cus - coarse_cells_y;

                    for(int fy = 0; fy < coarse_side_multiplier; fy++) {
                        for(int fx = 0; fx < coarse_side_multiplier; fx++) {
                            int y = coarse_side_multiplier * cy + fy;
                            int x = coarse_side_multiplier * cx + fx;
                            if(y == 0 || x == 0 || y == s->cells[1] - 1 || x == s->cells[0] - 1) {
                                continue;
                            }

                            int us = y * n + x;
                            int right = us + 1;
                            int left = us - 1;
                            int up = us + n;
                            int down = us - n;

                            coarse_solids[cus] += solids_buffer[us] / (fluid_type)pow(coarse_side_multiplier, 2);
                        }
                    }
                }
            }
            
        }
};

Scene* create_scene(int scene_type = STARTING_SCENARIO) {
    Scene* scene = new Scene();
    scene->setup_scene(scene_type);
    return scene;
}

struct LogRow {
    public:
        optional<int> iteration;
        optional<fluid_type> time;
        optional<fluid_type> avg_pressure;
        optional<fluid_type> min_pressure;
        optional<fluid_type> max_pressure;
        optional<fluid_type> avg_u_velocity;
        optional<fluid_type> avg_v_velocity;
        optional<fluid_type> avg_energy;
        
        bool ready() {
            return iteration.has_value() && time.has_value() && avg_pressure.has_value() && min_pressure.has_value() && max_pressure.has_value() && avg_u_velocity.has_value() && avg_energy.has_value();
        }

        void write(ofstream &log) {
            log << (iteration.has_value() ? to_string(iteration.value()) : "NULL");
            log << ",";
            log << (time.has_value() ? to_string(time.value()) : "NULL");
            log << ",";
            log << (avg_pressure.has_value() ? to_string(avg_pressure.value()) : "NULL");
            log << ",";
            log << (min_pressure.has_value() ? to_string(min_pressure.value()) : "NULL");
            log << ",";
            log << (max_pressure.has_value() ? to_string(max_pressure.value()) : "NULL");
            log << ",";
            log << (avg_u_velocity.has_value() ? to_string(avg_u_velocity.value()) : "NULL");
            log << ",";
            log << (avg_v_velocity.has_value() ? to_string(avg_v_velocity.value()) : "NULL");
            log << ",";
            log << (avg_energy.has_value() ? to_string(avg_energy.value()) : "NULL");
            log << "\n";
        }
};

/// adjust the vertical velocity of every cell according to a gravity function
void integrate(State* s, fluid_type dt) {
    fluid_type g_strength = s->g_strength;
    void (*g_pointer)(State*, fluid_type, fluid_type, int, int) = s->g_pointer;

    if(g_strength == 0.0) {
        return;
    }

    int n = s->cells[0];

    for(int y = 1; y < s->cells[1]-1; y++) {
        for(int x = 1; x < s->cells[0]-1; x++) {
            if(s->solids[x + n * y] != 0.0 && s->solids[x + n * (y - 1)] != 0.0) {
                //s->v_current[i*n + j] += g * dt;
                g_pointer(s, g_strength, dt, x, y);//(State* s, fluid_type strength, fluid_type dt, int x, int y)
            }
        }
    }
}


struct AABB {
    int x0;
    int y0;
    int x1;
    int y1;
};
///gauss-seidel iteration on an axis aligned bounding box
void solve_AABB_buffers(AABB aabb, unsigned int full_width, int* test_buffer, fluid_type* u_buffer, fluid_type* v_buffer, fluid_type* u_residuals, fluid_type* v_residuals, fluid_type* pressure, fluid_type* solids_buffer, fluid_type over_relax, fluid_type viscosity, fluid_type cp) {

    int n = full_width;

    for(int y = aabb.y0; y <= aabb.y1; y++) {
        for(int x = aabb.x0; x <= aabb.x1; x++) {
            if(test_buffer != nullptr) {
                test_buffer[y * n + x] = 0;
            }

            if(solids_buffer[y * n + x] == S_FIELD_SOLID) {
                continue;
            }

            unsigned int us = y * n + x;
            unsigned int left = us - 1;
            unsigned int right = us + 1;
            unsigned int up = (y + 1) * n + x;
            unsigned int down = (y - 1) * n + x;

            fluid_type sol = solids_buffer[us];
            fluid_type sx0 = solids_buffer[left];
            fluid_type sx1 = solids_buffer[right];
            fluid_type sy0 = solids_buffer[down];
            fluid_type sy1 = solids_buffer[up];

            sol = sx0 + sx1 + sy0 + sy1;

            if(sol == S_FIELD_SOLID) {
                continue;
            }

            fluid_type u_r0 = u_residuals == nullptr ? 0.0 : u_residuals[us];
            fluid_type v_r0 = v_residuals == nullptr ? 0.0 : v_residuals[us];
            fluid_type u_r1 = u_residuals == nullptr ? 0.0 : u_residuals[right];
            fluid_type v_r1 = v_residuals == nullptr ? 0.0 : v_residuals[up];

            fluid_type residual_divergence = u_r1 - u_r0 + v_r1 - v_r0;

            //every iteration the velocity in each edge decreases in magnitude by an amount proportional to 
            // the total divergence of both cells its attached to
            //Ax = (b = 0)
            fluid_type divergence = u_buffer[right] - u_buffer[us] + v_buffer[up] - v_buffer[us];
            fluid_type p = (residual_divergence - divergence) / sol;

            p *= over_relax;
            if(pressure != nullptr) {
                pressure[us] += cp * p;
            }

            u_buffer[us] += sx0 * (0 - p);
            u_buffer[right] -= sx1 * (0 - p);
            v_buffer[us] += sy0 * (0 - p);
            v_buffer[up] -= sy1 * (0 - p);
            if(isnan(u_buffer[us]) || isnan(u_buffer[right]) || isnan(v_buffer[us]) || isnan(v_buffer[up])) {
                printf("\n\tsolve_AABB_buffers just set a nan! x %i y %i u us %g u right %g v us %g v up %g u_r %g v_r %g sol %g", x, y, u_buffer[us], u_buffer[right], v_buffer[us], v_buffer[up], u_r0, v_r0, sol);
                paused = true;
                abort();
            }
        }
    }
}

// fine -> coarse
void restriction(State* s, int coarse_cells_x, int coarse_cells_y, int coarse_side_multiplier, fluid_type* u_buffer, fluid_type* v_buffer, fluid_type* u_residuals, fluid_type* v_residuals, fluid_type* solids_buffer, fluid_type* coarse_solids) {
    int n = s->cells[0];
    int cn = coarse_cells_y;
    for(int cy = 0; cy < coarse_cells_y; cy++) {
        for(int cx = 0; cx < coarse_cells_x; cx++) {
            fluid_type u_residual = 0.0; //the sum of the divergence of each of the cells we're downsampled from
            fluid_type v_residual = 0.0;
            int cus = cx + coarse_cells_y * cy;
            int cright = cus + 1;
            int cleft = cus - 1;
            int cup = cus + coarse_cells_y;
            int cdown = cus - coarse_cells_y;

            for(int fy = 0; fy < coarse_side_multiplier; fy++) {
                for(int fx = 0; fx < coarse_side_multiplier; fx++) {
                    int y = coarse_side_multiplier * cy + fy;
                    int x = coarse_side_multiplier * cx + fx;
                    if(y == 0 || x == 0 || y == s->cells[1] - 1 || x == s->cells[0] - 1) {
                        continue;
                    }

                    int us = y * n + x;
                    int right = us + 1;
                    int left = us - 1;
                    int up = us + n;
                    int down = us - n;

                    //coarse_solids[cx + coarse_cells_y * cy] += solids_buffer[us] / (fluid_type)pow(coarse_side_multiplier, 2);

                    //if(isnan(coarse_solids[cx + coarse_cells_y * cy])) {
                    //    printf("\n\trestriction had a nan in coarse_solids already!: cx %i cy %i fx %i fy %i x %i y %i solids[us] %g, pow(coarse_side_mult, 2) %g", 
                    //    cx,cy,fx,fy,x,y,solids_buffer[us],(fluid_type)pow(coarse_side_multiplier, 2));
                    //    abort();
                    //}
                    //coarse_solids[cx + coarse_cells_y * cy] += solids_buffer[us] / (fluid_type)pow(coarse_side_multiplier, 2);
                    //if(isnan(coarse_solids[cx + coarse_cells_y * cy])) {
                    //    printf("\n\trestriction just set a nan in coarse_solids: cx %i cy %i fx %i fy %i x %i y %i solids[us] %g, pow(coarse_side_mult, 2) %g", 
                    //    cx,cy,fx,fy,x,y,solids_buffer[us],(fluid_type)pow(coarse_side_multiplier, 2));
                    //    abort();
                    //}

                    if(solids_buffer[us] == S_FIELD_SOLID) {
                        continue;
                    }
                    fluid_type sol = solids_buffer[us];
                    fluid_type sx0 = solids_buffer[left];
                    fluid_type sx1 = solids_buffer[right];
                    fluid_type sy0 = solids_buffer[down];
                    fluid_type sy1 = solids_buffer[up];

                    sol = sx0 + sx1 + sy0 + sy1;

                    //Ax = b -> r = b - Ax
                    fluid_type divergence = u_buffer[right] - u_buffer[us] + v_buffer[up] - v_buffer[us];
                    fluid_type p = divergence / sol;

                    //u_residual -= sx0 * p;
                    //v_residual -= sy0 * p;
                    u_residuals[cx + coarse_cells_y * cy] -= sx0 * p;
                    if(cx + 1 + coarse_cells_y * cy < coarse_cells_x * coarse_cells_y) {
                        u_residuals[cx + 1 + coarse_cells_y * cy] += sx1 * p;
                    } 
                    v_residuals[cx + coarse_cells_y * cy] -= sy0 * p;
                    if(cx + coarse_cells_y * (cy + 1) < coarse_cells_x * coarse_cells_y) {
                        v_residuals[cx + coarse_cells_y * (cy + 1)] += sy1 * p;
                    }
                    //if(isnan(u_residuals[cx + coarse_cells_y * cy]) || isnan(u_residuals[cx + 1 + coarse_cells_y * cy]) || 
                    //    isnan(v_residuals[cx + coarse_cells_y & cy] || isnan(v_residuals[cx + coarse_cells_y * (cy + 1)]))) {
                    //    printf("\n\trestriction just set a nan! cx %i cy %i u residuals[us] %g u residuals[right] %g v residuals[us] %g v residuals[up]", cx, cy, 
                    //        u_residuals[cx + coarse_cells_y * cy], u_residuals[cx + 1 + coarse_cells_y * cy], v_residuals[cx + coarse_cells_y & cy], v_residuals[cx + coarse_cells_y * (cy + 1)]);
                    //    paused = true;
                    //    abort();
                    //}
                }
            }

            //u_residuals[cx + coarse_cells_y * cy] = u_residual;
            //v_residuals[cx + coarse_cells_y * cy] = v_residual;
            //if(isnan(u_residual) || isnan(v_residual)) {
            //    printf("\n\trestriction just set a nan! cx %i cy %i u residual %g v residual %g", cx, cy, u_residual, v_residual);
            //    paused = true;
            //}
        }
    }
}

//coarse -> fine via interpolation of the delta u and v solved using the coarse grid, 
// added onto the total u and v vels we iterated on before restriction
void prolongation(State* s, int coarse_cells_x, int coarse_cells_y, int coarse_side_multiplier, fluid_type* u_buffer, fluid_type* v_buffer, fluid_type* u_error, fluid_type* v_error, fluid_type* pressure_buffer, fluid_type cp, fluid_type* solids_buffer, fluid_type* coarse_solids) {
    int n = s->cells[0];
    int coarse_len = coarse_cells_x * coarse_cells_y;
    for(int y = 1; y < s->cells[1] - 1; y++) {
        for(int x = 1; x < s->cells[0] - 1; x++) {
            int us = x + y * n;
            int left = us - 1;
            int right = us + 1;
            int up = (y + 1) * n + x;
            int down = (y - 1) * n + x;
            if(solids_buffer[us] == S_FIELD_SOLID) {
                continue;
            }


            fluid_type raw_cy = (fluid_type)y / (fluid_type)coarse_side_multiplier;
            fluid_type raw_cx = (fluid_type)x / (fluid_type)coarse_side_multiplier;

            int cy = (int)floor(raw_cy);
            int cx = (int)floor(raw_cx);
            int cy_ahead = (int)ceil(raw_cy);
            int cx_ahead = (int)ceil(raw_cx);

            int cus = cx + coarse_cells_y * cy;
            int cright = cus + 1;
            int cleft = cus - 1;
            int cup = cus + coarse_cells_y;
            int cdown = cus - coarse_cells_y;

            //printf("\n\t\t--------PROLONG BEFORE C WEIGHTS");
            fluid_type cs_down = cy > 0 ? coarse_solids[cdown] : 1.0;
            fluid_type cs_up = cy < coarse_cells_y - 1 ? coarse_solids[cup] : 1.0;
            fluid_type cs_left = cx > 0 ? coarse_solids[cleft] : 1.0;
            fluid_type cs_right = cx < coarse_cells_x - 1 ? coarse_solids[cright] : 1.0;

            fluid_type cy_behind_weight = (cy_ahead - raw_cy) * solids_buffer[down];
            fluid_type cy_ahead_weight = (raw_cy - cy) * solids_buffer[up];
            
            fluid_type cx_behind_weight = (cx_ahead - raw_cx) * solids_buffer[left];
            fluid_type cx_ahead_weight = (raw_cx - cx) * solids_buffer[right];

            //printf("\n\tprolong x %i y %i cx %i cy %i raw cx %g raw cy %g cy behind %g ahead %g cx behind %g ahead %g",
            //    x,y,cx,cy,raw_cx,raw_cy, cy_behind_weight, cy_ahead_weight, cx_behind_weight, cx_ahead_weight);

            //printf("\n\t\t--------PROLONG BEFORE BUFFER FILLS");

            //pressure_buffer[us] += u_error[cus] + v_error[cus];
            fluid_type u_correction = (cx_behind_weight * u_error[cus] + cx_ahead_weight * (cright < coarse_cells_x ? u_error[cright] : u_error[cus]));
            fluid_type v_correction = (cy_behind_weight * v_error[cus] + cy_ahead_weight * (cup < coarse_cells_y ? v_error[cup] : v_error[cus]));
            //fluid_type divergence0 = u_buffer[right] - u_buffer[us] + v_buffer[up] - v_buffer[us];
            //fluid_type divergence1 = u_buffer[right] - u_buffer[us] + v_buffer[up] - v_buffer[us];
            
            //printf("\n\tprolong x %i y %i u %g -> %g, v %g -> %g", x,y,u_buffer[us], u_buffer[us] + u_correction, v_buffer[us], v_buffer[us] + v_correction);
            u_buffer[us] += u_correction;
            v_buffer[us] += v_correction;
            //if(isnan(u_buffer[us]) || isnan(v_buffer[us])) {
            //    printf("\n\tprolongation just set a nan! x %i y %i u us %g v us %g u error %g v error %g", x, y, u_buffer[us], v_buffer[us], u_error[us], v_error[us]);
            //    paused = true;
            //}
        }
    }
}

void fit_AABBs(AABB* AABBs, int rects, int rect_levels, int width, int height) {

    int rects_per_side = 2 * rect_levels;

    int rect_width = width / rects_per_side;
    int rect_height = height / rects_per_side;

    int* x_limits = new int[2 * rects_per_side];
    int* y_limits = new int[2 * rects_per_side];

    int px = 1;
    int py = 1;

    int x_rect_i = 0;
    int y_rect_i = 0;


    for(x_rect_i; x_rect_i < rects_per_side; x_rect_i++) {
        x_limits[2 * x_rect_i] = px;
        px += rect_width;
        x_limits[2 * x_rect_i + 1] = px - 1;
    }
    x_limits[2 * (rects_per_side - 1) + 1] = width - 2;

    for(int ry = 0; ry < rects_per_side; ry++) {
        y_limits[2 * ry] = py;
        py += rect_height;
        y_limits[2 * ry + 1] = py - 1;
    }
    y_limits[2 * (rects_per_side - 1) + 1] = height - 2;

    for(int ry = 0; ry < rects_per_side; ry++) {
        for(int rx = 0; rx < rects_per_side; rx++) {
            AABBs[rx + ry * rects_per_side].x0 = x_limits[2 * rx];
            AABBs[rx + ry * rects_per_side].x1 = x_limits[2 * rx + 1];

            AABBs[rx + ry * rects_per_side].y0 = y_limits[2 * ry];
            AABBs[rx + ry * rects_per_side].y1 = y_limits[2 * ry + 1];

            //printf("\n\t(%i, %i): (%i-%i, %i-%i)", rx, ry, AABBs[rx + ry * rects_per_side].x0, AABBs[rx + ry * rects_per_side].x1, AABBs[rx + ry * rects_per_side].y0, AABBs[rx + ry * rects_per_side].y1);
        }
    }

    delete [] x_limits;
    delete [] y_limits;
}

fluid_type calc_avg_divergence(State* s) {
    int n = s->cells[0];
    fluid_type* u_buffer = s->u_current;
    fluid_type* v_buffer = s->v_current;
    fluid_type* solids_buffer = s->solids;
    fluid_type avg_divergence = 0.0;
    for(int y = 1; y < s->cells[1] - 1; y++) {
        for(int x = 1; x < s->cells[0] - 1; x++) {

            unsigned int us = y * n + x;
            unsigned int left = us - 1;
            unsigned int right = us + 1;
            unsigned int up = (y + 1) * n + x;
            unsigned int down = (y - 1) * n + x;

            fluid_type sol = solids_buffer[us];
            fluid_type sx0 = solids_buffer[left];
            fluid_type sx1 = solids_buffer[right];
            fluid_type sy0 = solids_buffer[down];
            fluid_type sy1 = solids_buffer[up];

            sol = sx0 + sx1 + sy0 + sy1;

            if(sol == S_FIELD_SOLID) {
                continue;
            }

            avg_divergence -= (sx1*u_buffer[right] - sx0*u_buffer[us] + sy1*v_buffer[up] - sy0*v_buffer[us]) / sol;
        }
    }
    avg_divergence /= s->num_cells;
    return abs(avg_divergence);
}

void solve_incompressability_multithread_multigrid(State* s, int num_iterations, fluid_type dt) {
    int n = s->cells[0];
    fluid_type cp = s->density * s->h / dt;

    
    int rect_levels = 2; //0 = no recursion, 1 = world rect is split into 4 subrects, 2 = 1 and each subrect is split into 4 subsubrects. every rect gets a thread
    int rects = (int)pow(4, rect_levels);
    AABB* AABBs = new AABB[rects];
    //vector<thread> threads;
    thread** threads = new thread*[rects];

    fit_AABBs(AABBs, rects, rect_levels, s->cells[0], s->cells[1]);

    if(STARTING_COARSE_GRID_BUDGET_DIVISOR == 0 || STARTING_COARSE_GRID_SIDE_LENGTH_POWER == 0) {

        unsigned int full_width = s->cells[0];
        int* test_buffer = s->debug_touch_buffer;

        for(int i = 0; i < s->num_cells; i++) {
            test_buffer[i] = 1;
        }

        printf("\n-------BEFORE FIRST FINE ITER");
        fluid_type avg_divergence = calc_avg_divergence(s);
        printf("\n\tbefore iter, avg divergence is %g", avg_divergence);

        //first fine iteration portion
        for(int iter = 0; iter < num_iterations; iter++) {
            
            for(int i = 0; i < rects; i++) {
                threads[i] = new thread(solve_AABB_buffers, AABBs[i], full_width, test_buffer, s->u_current, s->v_current, nullptr, nullptr, s->pressure, s->solids, s->over_relaxation, s->viscosity, cp);
            }

            for(int i = 0; i < rects; i++) {
                threads[i]->join();
                delete threads[i];
            }
        }
        fluid_type avg_divergence2 = calc_avg_divergence(s);
        printf("\n\tafter iter, avg divergence is %g (%g)", avg_divergence2, avg_divergence2 - avg_divergence);
        delete [] AABBs;
        delete [] threads;
        return;
    }
   
    //we split our iterations into 3 regions: the first fine iterations, the coarse iterations, and the second fine iterations
    //(1/n) of num_iterations is in the coarse region, the remaining is evenly split between both fine regions
    //we spend the first fine iterations on fine iterations to reduce the error a little bit,
    //then we downsample the u and v velocities into a grid where the cells are 2^m times larger (and thus theres 2^-m as many cells to iterate through)
    //then we do coarse iterations on that grid, and when we're done we interpolate each coarse cell into the appropiate number of 
    //fine cells again, and then we do the remaining iteration budget on the fine grid again

    //the coarse grid will take up 1/this of the total iteration budget
    fluid_type coarse_grid_iteration_budget_divisor = STARTING_COARSE_GRID_BUDGET_DIVISOR;
    //how many iterations we do on the coarse grid
    int coarse_iterations = coarse_grid_iteration_budget_divisor > 0.0 ? (int)floor((fluid_type)num_iterations / coarse_grid_iteration_budget_divisor) : 0;
    //how many iterations we do on the fine grid in total
    int fine_iterations = num_iterations - coarse_iterations;
    int first_fine_iterations = (int)(fine_iterations / 2);
    int second_fine_iterations = fine_iterations - first_fine_iterations;

    int coarse_grid_side_length_power = STARTING_COARSE_GRID_SIDE_LENGTH_POWER; //coarse cells are 2^this times as big per side, so 2^-this as many of them
    int coarse_grid_side_length_mult = (int)pow(2, coarse_grid_side_length_power);

    int coarse_grid_rect_levels_modifier = 1;
    int coarse_grid_rect_levels = rect_levels - coarse_grid_rect_levels_modifier;
    int coarse_rects = (int)pow(4, coarse_grid_rect_levels);
    AABB* coarse_AABBs = new AABB[coarse_rects];
    thread** coarse_threads = new thread*[coarse_rects];

    int coarse_cells_x = s->cells[0] / coarse_grid_side_length_mult;
    int coarse_cells_y = s->cells[1] / coarse_grid_side_length_mult;
    int coarse_cells = coarse_cells_x * coarse_cells_y;
    int fine_cells_per_coarse = (int)pow(coarse_grid_side_length_mult, 2);
    fluid_type coarse_h = s->h * coarse_grid_side_length_mult;

    fit_AABBs(coarse_AABBs, coarse_rects, coarse_grid_rect_levels, coarse_cells_x, coarse_cells_y);

    // error terms. starts at 0, then solved in the coarse grid iteratively, then 
    // in prolongation we add the interpolated values of these at each coarse cell into the
    // actual u and v velocities. THEN we do a final gauss seidel fine iteration round
    // to get the final approximated 0 divergence velocity field
    fluid_type* u_coarse = s->u_coarse;
    fluid_type* v_coarse = s->v_coarse;
    // residual terms
    fluid_type* u_residuals = s->u_residuals;
    fluid_type* v_residuals = s->v_residuals;
    //A matrix but restricted into the coarse grid
    fluid_type* solids_coarse = s->solids_coarse;

    for(int ci = 0; ci < s->coarse_cells; ci++) {
        u_coarse[ci] /= 2; //to keep them from growing out of control
        v_coarse[ci] /= 2;
    }

    //i dont think i need this
    //fluid_type* pressure_coarse = new fluid_type[coarse_cells];

    unsigned int full_width = s->cells[0];
    unsigned int width = s->cells[0] / 2;
    unsigned int height = s->cells[1] / 2;

    unsigned int width2 = (width / 2);
    unsigned int height2 = height / 2;

    int* test_buffer = s->debug_touch_buffer;

    for(int i = 0; i < s->num_cells; i++) {
        test_buffer[i] = 1;
    }

    //printf("\n-------BEFORE FIRST FINE ITER");

    //fluid_type avg_divergence = calc_avg_divergence(s);
    //printf("\n\tbefore 1st fine iter, avg divergence is %g", avg_divergence);
    //first fine iteration portion
    for(int iter = 0; iter < first_fine_iterations; iter++) {
        
        for(int i = 0; i < rects; i++) {
            threads[i] = new thread(solve_AABB_buffers, AABBs[i], full_width, test_buffer, s->u_current, s->v_current, nullptr, nullptr, s->pressure, s->solids, s->over_relaxation, s->viscosity, cp);
        }

        for(int i = 0; i < rects; i++) {
            threads[i]->join();
            delete threads[i];
        }
    }

    //fluid_type avg_divergence2 = calc_avg_divergence(s);
    //printf("\n\tafter 1st fine iter, avg divergence is %g (%g)", avg_divergence2, avg_divergence2 - avg_divergence);
    //avg_divergence = avg_divergence2;

    //printf("\n-------BEFORE RESTRICTION");
    restriction(s, coarse_cells_x, coarse_cells_y, coarse_grid_side_length_mult, s->u_current, s->v_current, u_coarse, v_coarse, s->solids, solids_coarse);


    //printf("\n-------BEFORE COARSE ITER");
    //coarse iterations
    AABB coarse_AABB;
    coarse_AABB.x0 = 1;
    coarse_AABB.x1 = coarse_cells_x - 2;
    coarse_AABB.y0 = 1;
    coarse_AABB.y1 = coarse_cells_y - 2;
    for(int iter = 0; iter < coarse_iterations; iter++) {
        //solve_AABB_buffers(coarse_AABB, coarse_cells_x, nullptr, u_coarse, v_coarse, u_residuals, v_residuals, nullptr, solids_coarse, s->over_relaxation, s->viscosity, cp * coarse_h / s->h);
        for(int i = 0; i < coarse_rects; i++) {
            coarse_threads[i] = new thread(solve_AABB_buffers, coarse_AABBs[i], coarse_cells_x, nullptr, u_coarse, v_coarse, u_residuals, v_residuals, nullptr, solids_coarse, 1.5, s->viscosity, cp * coarse_h / s->h);
        }
        
        for(int i = 0; i < coarse_rects; i++) {
            coarse_threads[i]->join();
            delete coarse_threads[i];
        }
    }

    //printf("\n-------BEFORE PROLONGATION");
    prolongation(s, coarse_cells_x, coarse_cells_y, coarse_grid_side_length_mult, s->u_current, s->v_current, u_coarse, v_coarse, s->pressure, cp, s->solids, solids_coarse);

    //avg_divergence2 = calc_avg_divergence(s);
    //printf("\n\tafter coarse iter and prolongation, avg divergence is %g (%g)", avg_divergence2, avg_divergence2 - avg_divergence);
    //avg_divergence = avg_divergence2;
    //now we need to set u and v and pressure in the fine grid to an interpolation of the averages of the coarse grid values

    //printf("\n-------BEFORE SECOND FINE ITER");
    //second fine iterations
    for(int iter = 0; iter < second_fine_iterations; iter++) {
        
        for(int i = 0; i < rects; i++) {
            threads[i] = new thread(solve_AABB_buffers, AABBs[i], full_width, test_buffer, s->u_current, s->v_current, nullptr, nullptr, s->pressure, s->solids, s->over_relaxation, s->viscosity, cp);
        }

        for(int i = 0; i < rects; i++) {
            threads[i]->join();
            delete threads[i];
        }
    }

    //avg_divergence2 = calc_avg_divergence(s);
    //printf("\n\tafter 2nd fine iter, avg divergence is %g (%g)", avg_divergence2, avg_divergence2 - avg_divergence);

    

    delete [] AABBs;
    delete [] threads;
    delete [] coarse_AABBs;
    delete [] coarse_threads;
}

void solve_incompressability_multithread(State* s, int num_iterations, fluid_type dt) {
    int n = s->cells[0];
    fluid_type cp = s->density * s->h / dt;

    
    int rect_levels = 2; //0 = no recursion, 1 = world rect is split into 4 subrects, 2 = 1 and each subrect is split into 4 subsubrects. every rect gets a thread
    int rects = (int)pow(4, rect_levels);
    AABB* AABBs = new AABB[rects];
    //vector<thread> threads;
    thread** threads = new thread*[rects];
    vector<AABB> AABBs2;
    //for(int i = 0; i < rects; i++) {
    //    AABBs2.push_back(AABBs[i]);
    //}

    fit_AABBs(AABBs, rects, rect_levels, s->cells[0], s->cells[1]);

    unsigned int full_width = s->cells[0];
    unsigned int width = s->cells[0] / 2;
    unsigned int height = s->cells[1] / 2;

    unsigned int width2 = (width / 2);
    unsigned int height2 = height / 2;

    int* test_buffer = s->debug_touch_buffer;

    for(int i = 0; i < s->num_cells; i++) {
        test_buffer[i] = 1;
    }

    for(int iter = 0; iter < num_iterations; iter++) {

        //solve_AABB_buffers(AABB aabb, unsigned int full_width, int* test_buffer, fluid_type* u_buffer, fluid_type* v_buffer, fluid_type* pressure, fluid_type* solids_buffer, fluid_type over_relax, fluid_type viscosity, fluid_type cp) {
        //for(AABB &aabb : AABBs2) {
        //    threads.emplace_back(std::ref(solve_AABB_buffers), aabb, full_width, test_buffer, s->u_current, s->v_current, s->pressure, s->solids, s->over_relaxation, cp);
        //}
        
        for(int i = 0; i < rects; i++) {
            //threads[i] = new thread(solve_AABB_buffers, AABBs[i], full_width, test_buffer, s->u_current, s->v_current, s->pressure, s->solids, s->over_relaxation, s->viscosity, cp);
        }

        for(int i = 0; i < rects; i++) {
            threads[i]->join();
            delete threads[i];
        }

    }
    
    if(false) {

        for(int y = 1; y < s->cells[1] - 1; y++) {
            for(int x = 1; x < s->cells[0] - 1; x++) {
                if(test_buffer[y * s->cells[0] + x] == 1) {
                    string x_offset = "";
                    string y_offset = "";
                    //if(x == )
                    printf("\n\t%i %i untouched", x, y);
                }
            }
        }

    }

    delete [] AABBs;

}

void solve_incompressability(State* s, int num_iterations, fluid_type dt) {
    int n = s->cells[0];
    fluid_type cp = s->density * s->h / dt;

    for(int iter = 0; iter < num_iterations; iter++) {
        for(int y = 1; y < s->cells[1]-1; y++) {
            for(int x = 1; x < s->cells[0]-1; x++) {

                if(s->solids[y * n + x] == S_FIELD_SOLID) {
                    continue;
                }

                unsigned int us = y * n + x;
                unsigned int left = us - 1;
                unsigned int right = us + 1;
                unsigned int up = (y + 1) * n + x;
                unsigned int down = (y - 1) * n + x;

                fluid_type sol = s->solids[us];
                fluid_type sx0 = s->solids[left];
                fluid_type sx1 = s->solids[right];
                
                fluid_type sy0 = s->solids[down];
                fluid_type sy1 = s->solids[up];
                sol = sx0+sx1+sy0+sy1;
                if(sol == S_FIELD_SOLID) {
                    continue;
                }

                fluid_type divergence = s->u_current[right] - s->u_current[us] + s->v_current[up] - s->v_current[us];
                fluid_type p = -divergence / sol;

                p *= s->over_relaxation;
                s->pressure[us] += cp * p;


                s->u_current[us] -= sx0 * p / s->viscosity;
                s->u_current[right] += sx1 * p / s->viscosity;
                s->v_current[us] -= sy0 * p / s->viscosity;
                s->v_current[up] += sy1 * p / s->viscosity;
            }
        }
    }
}

/// @brief sets boundary velocities equal to the adjacent cell closest to the center's velocity
/// @param s 
void extrapolate(State* s) {
    int n = s->cells[0];
    for(int x = 0; x < s->cells[0]; x++) {
        s->u_current[x] = s->u_current[n + x];
    }
    for(int x = 0; x < s->cells[0]; x++) {
        s->u_current[(s->cells[1] - 1) * n + x] = s->u_current[(s->cells[1] - 2) * n + x];
    }

    for(int y = 0; y < s->cells[1]; y++) {
        s->v_current[y * n] = s->v_current[y * n + 1];
        s->v_current[y * n + s->cells[0] - 1] = s->v_current[y * n + s->cells[0] - 2];
    }
}

/**
 * x and y are distance normalized x and y coordinates - dt * velocities at that point.
 * x and y are a point in space, which lies between 2-4 points corresponding to the actual indices of that field, which we want to average
 * to get the returned quantity at the given point. what points the given point is between depends on the field: u velocities are stored on the
 * horizontal edges between cells, so we want to take a weighted average between u_current[closest index to our left] and u_current[that index + 1]
 * while for v velocities we need to lerp between v_current[closest index below us] and v_current[that index + ]
 */
fluid_type sample_field(State* s, fluid_type x, fluid_type y, int field) {
    int n = s->cells[0];
    fluid_type h = s->h;
    fluid_type h1 = 1.0 / h;
    fluid_type h2 = 0.5 * h;

    x = max(min(x, s->cells[0] * h), h);
    y = max(min(y, s->cells[1] * h), h);

    fluid_type dx = 0.0;
    fluid_type dy = 0.0;

    fluid_type* f;

    switch(field) {
        case U_FIELD:
            f = s->u_current;
            dy = h2;
            break;
        case V_FIELD:
            f = s->v_current;
            dx = h2;
            break;
        case S_FIELD:
            f = s->m_current;
            dx = h2;
            dy = h2;
            break;
        case T_FIELD:
            f = s->temperature_current;
            dx = h2;
            dy = h2;
            break;
    }

    //closest edge on our left
    int x0 = max(min((int)floor((x - dx) / h), s->cells[0] - 1), 0); // (x - change in x) / h
    fluid_type tx = ((x-dx) - x0*h) / h;
    //closest edge on our right
    int x1 = max(min(x0 + 1, s->cells[0]-1), 0);

    //closest edge below us
    int y0 = max(min((int)floor((y-dy) / h), s->cells[1]-1), 0);
    fluid_type ty = ((y-dy) - y0*h) / h;
    //closest edge above us
    int y1 = max(min(y0 + 1, s->cells[1]-1), 0);

    fluid_type sx = 1.0 - tx;
    fluid_type sy = 1.0 - ty;

    //printf("\n\tsample_field x %g y %g x0 %i x1 %i y0 %i y1 %i", x, y, x0, x1, y0, y1);

    fluid_type val = sx*sy * f[x0+n*y0] + tx*sy * f[x1+n*y0] + tx*ty * f[x1+n*y1] + sx*ty * f[x0 + n * y1];

    return val;
}

/*

			var x0 = Math.min(Math.floor((x-dx)*h1), this.numX-1);
			var tx = ((x-dx) - x0*h) * h1;
			var x1 = Math.min(x0 + 1, this.numX-1);
			
			var y0 = Math.min(Math.floor((y-dy)*h1), this.numY-1);
			var ty = ((y-dy) - y0*h) * h1;
			var y1 = Math.min(y0 + 1, this.numY-1);

			var sx = 1.0 - tx;
			var sy = 1.0 - ty;

			var val = sx*sy * f[x0*n + y0] +
				tx*sy * f[x1*n + y0] +
				tx*ty * f[x1*n + y1] +
				sx*ty * f[x0*n + y1];
			
			return val;
*/

/// u = horizontal velocities at the vertical (down to up) edges of each cell
/// u[index of cell] = horizontal velocity at the leftmost edge of that cell
/// this gives the average of the horizontal velocities of our left and right edges and the left and right edges of the cell below us
fluid_type avg_u(State* s, int i, int j) {
    int n = s->cells[0];
    fluid_type u = 0.25 * (s->u_current[i+ n*(j-1)] + s->u_current[i + n * j] + s->u_current[(i+1) + n*(j-1)] + s->u_current[(i-1)+ n*(j+1)]);
    //fluid_type u = ((1.0/total_weight) * s->u_current[j * n + i] + (coeff / total_weight) * (s->u_current[(j-1) * n + i] + s->u_current[(j-1) * n + i + 1] + s->u_current[j * n + i + 1]));
    return u;
}

/// v = vertical velocities at the horizontal (left to right) edges of each cell
/// v[index of cell] = vertical velocity at the lower edge of that cell.
/// this gives the average of the vertical velocities of our upper and lower edges and the upper and lower edges of the cell to our left
fluid_type avg_v(State* s, int i, int j) {
    int n = s->cells[0];
    fluid_type v = 0.25 * (s->v_current[(i-1) +n*j] + s->v_current[i + n * j] + s->v_current[(i-1) + n*(j+1)] + s->v_current[i + n*(j+1)]);
    //fluid_type v = (1.0 / total_weight) * (s->v_current[j * n + i - 1]) + (coeff / total_weight) * (s->v_current[j*n + i] + s->v_current[(j+1) * n + i - 1] + s->v_current[(j+1)*n + i]);
    return v;

}

void advect_vel(State* s, fluid_type dt) {


    int n = s->cells[0];
    fluid_type h = s->h;
    fluid_type h2 = 0.5 * h;

    fluid_type kinematic_viscosity = s->viscosity / s->density; //idk why i use this, we want mu which is the dynamic viscosity

    for(int i = 0; i < s->num_cells; i++) {//why do i do this this is dumb
        s->u_next[i] = s->u_current[i];
        s->v_next[i] = s->v_current[i];
    }

    fluid_type* sol = s->solids;

    fluid_type* u_next = s->u_next;
    fluid_type* v_next = s->v_next;
    fluid_type* u_current = s->u_current;
    fluid_type* v_current = s->v_current;
    fluid_type viscosity = s->viscosity;

    for(int y = 1; y < s->cells[1] - 1; y++) {
        for(int x = 1; x < s->cells[0] - 1; x++) { 
            int us = x + n * y;
            if(sol[us] == S_FIELD_SOLID) {
                continue;
            }
            int left = us - 1;
            int up = x + n * (y+1);
            int right = x + 1 + n * y;
            int down = x + n * (y - 1);

            //printf("\n%i, %i", x, y);
            fluid_type vleft = sol[left] != S_FIELD_SOLID ? s->v_current[left] : s->v_current[us];
            fluid_type vright = sol[right] != S_FIELD_SOLID ? s->v_current[right] : s->v_current[us];
            fluid_type vup = sol[up] != S_FIELD_SOLID ? s->v_current[up] : s->v_current[us];
            fluid_type vdown = sol[down] != S_FIELD_SOLID ? s->v_current[down] :  s->v_current[us];


            fluid_type uleft = sol[left] != S_FIELD_SOLID ? s->u_current[left] : s->u_current[us];
            fluid_type uright = sol[right] != S_FIELD_SOLID ? s->u_current[right] : s->u_current[us];
            fluid_type uup = sol[up] != S_FIELD_SOLID ? s->u_current[up] : s->u_current[us];
            fluid_type udown = sol[down] != S_FIELD_SOLID ? s->u_current[down] : s->u_current[us];

            //u component
            if(sol[left] != S_FIELD_SOLID) {
                fluid_type xi = x * h;
                fluid_type yi = y * h + h2;
                fluid_type u = s->u_current[us]; 
                fluid_type u_i = u;

                //fluid_type v = 0.2 * (vleft + s->v_current[us] + vup + vright + vdown); 
                //fluid_type v = 0.25 * (s->v_current[down] + s->v_current[us] + s->v_current[(y-1)*n+x+1] + s->v_current[y*n+x+1]);//avg_v(s, x, y);
                fluid_type v = avg_v(s, x, y);

                xi -= dt * u;
                yi -= dt * v;
                //printf("u sample");
                u = sample_field(s, xi, yi, U_FIELD);

                //printf("\n\t advect_vel u comp: x %i y %i u %g v %g u initial %g xi %g yi %g", x, y, u, v, u_i, xi, yi);
                s->u_next[us] = u;//(1.0 / (viscosity + 1.0)) * u + (viscosity / (viscosity + 1.0)) * 0.25 * (uleft + uright + uup + udown);
            }

            //v component
            if(sol[down] != S_FIELD_SOLID) {
                fluid_type xi = x * h + h2;
                fluid_type yi = y * h;

                //fluid_type u = 0.2 * (uleft + s->u_current[us] + uup + uright + udown); //0.25 * (s->u_current[left] + s->u_current[us] + s->u_current[(y+1)*n+x-1] + s->u_current[(y-1)*n+x+1]); //avg_u(s, x, y);
                fluid_type u = avg_u(s,x,y);
                fluid_type v = s->v_current[us];
                fluid_type v_i = v;

                xi -= dt * u;
                yi -= dt * v;

                //printf("v sample");
                v = sample_field(s, xi, yi, V_FIELD);
                //printf("\n\t advect_vel v comp: x %i y %i u %g v %g v initial %g xi %g yi %g", x, y, u, v, v_i, xi, yi);
                s->v_next[us] = v;//(1.0 / (viscosity + 1.0)) * v + (viscosity / (viscosity + 1.0) * (0.25 * (vleft + vright + vup + vdown)));//max(v, s->v_source[us]);
                //s->u_current[us] = (1.0 / (viscosity + 1.0)) * u_next[us] + (viscosity / (viscosity + 1.0)) * (0.25 * (uleft + uright + uup + udown));
                //s->v_current[us] = (1.0 / (viscosity + 1.0)) * v_next[us] + (viscosity / (viscosity + 1.0)) * (0.25 * (vleft + vright + vup + vdown));
            }


        }
    }

    //for(int i = 0; i < s->num_cells; i++) {//why do i do this this is dumb
    //    s->u_current[i] = s->u_next[i];
    //    s->v_current[i] = s->v_next[i];
    //}

    for(int y = 2; y < s->cells[1] - 1; y++) {
        for(int x = 2; x < s->cells[0] - 1; x++) {
            int us = x + y * n;
     
            if(sol[us] == S_FIELD_SOLID) {
                continue;
            }

            int left = us - 1;
            int right = us + 1;
            int up = x + n * (y+1);
            int down = x + n * (y-1);

            fluid_type uup = sol[up] != S_FIELD_SOLID ? u_next[up] : 0.0;
            fluid_type udown = sol[down] != S_FIELD_SOLID ? u_next[down] : 0.0;
            fluid_type uleft = sol[left] != S_FIELD_SOLID ? u_next[left] : 0.0;
            fluid_type uright = sol[right] != S_FIELD_SOLID ? u_next[right] : 0.0; 

            fluid_type vup = sol[up] != S_FIELD_SOLID ? v_next[up] : 0.0;
            fluid_type vdown = sol[down] != S_FIELD_SOLID ? v_next[down] : 0.0;
            fluid_type vleft = sol[left] != S_FIELD_SOLID ? v_next[left] : 0.0;
            fluid_type vright = sol[right] != S_FIELD_SOLID ? v_next[right] : 0.0; 
            //fluid_type x_stress = s->viscosity * (s->u_current[us] - s->u_current[x + (y-1) * n]) / s->h;
            //fluid_type y_stress = s->viscosity * (s->v_current[us] - s->v_current[left]) / s->h;

            s->u_current[us] = (1.0 / (viscosity + 1.0)) * u_next[us] + (viscosity / (viscosity * 1.05 + 1.0)) * (0.25 * (uleft + uright + uup + udown));
            s->v_current[us] = (1.0 / (viscosity + 1.0)) * v_next[us] + (viscosity / (viscosity * 1.05 + 1.0)) * (0.25 * (vleft + vright + vup + vdown));
        }

    }
}

void solve_viscid_flow(State* s, fluid_type dt) {
    int n = s->cells[0];
    fluid_type mu = s->viscosity;
    fluid_type* u_current = s->u_current;
    fluid_type* v_current = s->v_current;
    //having this as a separate step may lead to nonzero divergence but whatever fuck it
    for(int y = 1; y < s->cells[1]; y++) {
        for(int x = 1; x < s->cells[0]; x++) {
            int us = x + y * n;
            if(s->solids[us] == S_FIELD_SOLID) {
                continue;
            }

            int left = us - 1;
            int right = us + 1;
            int up = x + (y + 1) * n;
            int down = x + (y - 1) * n;

            //u is averaged between us and every neighbor
            fluid_type avg_u_x = u_current[us] - (1.0/2.0) * (u_current[left] + u_current[right]);
            fluid_type avg_u_y = u_current[us] - (1.0/2.0) * (u_current[up] + u_current[down]);

            u_current[us] += mu * (avg_u_x + avg_u_y);

            //v
            fluid_type avg_v_x = v_current[us] - (1.0/2.0) * (v_current[left] + v_current[right]);
            fluid_type avg_v_y = v_current[us] - (1.0/2.0) * (v_current[up] + v_current[down]);
            
            v_current[us] += mu * (avg_v_x + avg_v_y);
        }
    }
}

void advect_smoke(State* s, fluid_type dt) {
    for(int i = 0; i < s->num_cells; i++) {
        s->m_next[i] = s->m_current[i];
        //s->temperature_next[i] = s->temperature_current[i];
    }

    int n = s->cells[0];
    fluid_type h = s->h;
    fluid_type h2 = 0.5 * h;

    fluid_type diffusion_coeff = 0.3 * dt;

    double avg_smoke = 0.0;

    for(int y = 1; y < s->cells[1] - 1; y++) {
        for(int x = 1; x < s->cells[0] - 1; x++) {

            int us = x + n * y;
            int right = us + 1;
            int up = (x + n * (y + 1));
            if(s->solids[us] != S_FIELD_SOLID) {
                fluid_type u = (s->u_current[us] + s->u_current[right]) * 0.5;
                fluid_type v = (s->v_current[us] + s->v_current[up]) * 0.5;
                
                fluid_type xi = x * h + h2 - dt * u;
                fluid_type yi = y * h + h2 - dt * v;

                s->m_next[us] = sample_field(s, xi, yi, S_FIELD);
                avg_smoke += s->m_next[us];
                //s->temperature_next[us] = (sample_field(s, xi - diffusion_coeff, yi - diffusion_coeff, T_FIELD) + sample_field(s, xi + diffusion_coeff, yi + diffusion_coeff, T_FIELD)) / 2.0;
            }
        }
    }

    avg_smoke /= s->num_cells;
    //printf("\n\t\taverage smoke is %g", avg_smoke);

    for(int i = 0; i < s->num_cells; i++) {
        s->m_current[i] = s->m_next[i];
        //s->temperature_current[i] = s->temperature_next[i];
    }
}

void sim_step(State* s, int iterations, fluid_type time, fluid_type dt, LogRow* log_row, int scene_type, int workers_per_dim, unsigned int shader_program, unsigned int compute_program) {

    zero_out(s->debug_data_buffer, s->num_cells);
    zero_out(s->debug_touch_buffer, s->num_cells);
    zero_out(s->u_residuals, s->coarse_cells);
    zero_out(s->v_residuals, s->coarse_cells);

    //integration = add gravity acceleration to all velocities
    clock_t start = clock();
    integrate(s, dt);
    if(print_perf){printf("\n\t\tintegration took %g milliseconds", (fluid_type)(clock() - start) / (CLOCKS_PER_SEC / 1000.0));}

    if(s->clear_pressure == true) {
        for(int i = 0; i < s->num_cells; i++) {
            s->pressure[i] = 0.0;
        }
    }

    //make sure no fluid is created or destroyed
    start = clock();
    //glUseProgram(compute_program);
    //for(int i = 0; i < iterations; i++) {
    //    glDispatchCompute((GLuint)s->cells[0] / workers_per_dim, (GLuint)s->cells[1] / workers_per_dim, 1);
    //    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    //}
    //glUseProgram(shader_program);
    solve_incompressability_multithread_multigrid(s, s->num_divergence_iterations, dt);
    if(print_perf) {

        printf("\n\t\tsolve_incompressability() took %g milliseconds (%i iterations)", (fluid_type)(clock() - start) / (CLOCKS_PER_SEC / 1000.0), s->num_divergence_iterations);
    } 

    // fixes boundary conditions
    start = clock();
    extrapolate(s);
    if(print_perf){printf("\n\t\textrapolate() took %g milliseconds", (fluid_type)(clock() - start) / (CLOCKS_PER_SEC / 1000.0));}

    //samples the velocity of each point, goes to the point AWAY from that velocity
    //and sets the original points velocity to an average of the points at the far velocity
    start = clock();
    advect_vel(s, dt);
    if(print_perf){printf("\n\t\tadvect_vel() took %g milliseconds", (fluid_type)(clock() - start) / (CLOCKS_PER_SEC / 1000.0));}

    //move smoke like you did fluids
    start = clock();
    advect_smoke(s, dt);
    if(print_perf){printf("\n\t\tadvect_smoke() took %g milliseconds", (fluid_type)(clock() - start) / (CLOCKS_PER_SEC / 1000.0));}

    fluid_type min_pressure = 100000.0;
    fluid_type max_pressure = -1000000.0;
    fluid_type avg_pressure = 0.0;

    fluid_type avg_u_velocity = 0.0;
    fluid_type avg_v_velocity = 0.0;
    fluid_type avg_energy = 0.0;

    int n = s->cells[0];
    for(int y = 1; y < s->cells[1] - 1; y++) {
        for(int x = 1; x < s->cells[0] - 1; x++) { 
            int us = x + n * y;
            int right = us + 1;
            int left = us - 1;
            int up = (x + n * (y + 1));
            int down = x + n * (y - 1);

            fluid_type p_at = s->pressure[us];
            if(p_at > max_pressure) {
                max_pressure = p_at;
            }
            if(p_at < min_pressure) {
                min_pressure = p_at;
            }

            avg_pressure += p_at;

            fluid_type net_u = s->u_current[us] - s->u_current[left];
            fluid_type net_v = s->v_current[us] - s->v_current[down];
            avg_u_velocity += abs(net_u);
            avg_v_velocity += abs(net_v);
            
            fluid_type cell_mass = s->density * powf(s->h, 2.0);
            fluid_type energy = 0.5 * cell_mass * (powf(net_u, 2.0) + powf(net_v, 2.0));
            avg_energy += energy;
        }
    }

    avg_pressure = avg_pressure / s->num_cells;
    avg_u_velocity /= s->num_cells;
    avg_v_velocity /= s->num_cells;
    avg_energy /= s->num_cells;

    log_row->avg_pressure = avg_pressure;
    log_row->min_pressure = min_pressure;
    log_row->max_pressure = max_pressure;
    log_row->avg_u_velocity = avg_u_velocity;
    log_row->avg_v_velocity = avg_v_velocity;
    log_row->avg_energy = avg_energy;
    
    string clearing = "false";
    if(s->clear_pressure == true) {
        clearing = "true";
    }
    if(print_perf){printf("\n\t\tclear_pressure is %s, pressure from %g - %g, average is %g", clearing.c_str(), min_pressure, max_pressure, avg_pressure);}
    if(print_perf){printf("\n\t\taverage u vel is %g, avg v vel is %g, avg energy is %g", avg_u_velocity, avg_v_velocity, avg_energy);}

    if(print_perf) {
        printf("\n\t----------------------------------------------------------------------");
    }
}

bool is_sim_done(State* s, int iterations, fluid_type time, fluid_type max_time) {
    if(time > max_time) {
        return true;
    }
    return false;
}

void getSciColor(unsigned char* color, fluid_type val, fluid_type min_val, fluid_type max_val) {
    val = min(max(val, min_val), max_val - 0.00001f);
    fluid_type distance = max_val - min_val;
    if(distance == 0.0) {
        val = 0.5;
    } else {
        val = (val - min_val) / distance;
    }

    fluid_type m = 0.25;
    int num = (int)floor(val / m);
    fluid_type s = (val - num * m) / m;

    fluid_type r = 0.0;
    fluid_type g = 0.0; 
    fluid_type b = 0.0;

    switch(num) {
        case 0:
            r = 0.0; g = s; b = 1.0;
            break;
        case 1:
            r = 0.0; g = 1.0; b = 1.0 - s;
            break;
        case 2:
            r = s; g = 1.0; b = 0.0;
            break;
        case 3:
            r = 1.0; g = 1.0 - s; b = 0.0;
            break;
    }

    color[0] = 255 * r;
    color[1] = 255 * g;
    color[2] = 255 * b;
    color[3] = 255;
}


void error_callback(int error, const char* description) {
    fprintf(stderr, "Error: %s\n", description);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        printf("\nesc pressed, closing.");
        glfwSetWindowShouldClose(window, GLFW_TRUE);
        return;
    }

    if(action != GLFW_PRESS) {
        return;
    }

    if(change_mode_active == true) {
        if(key == GLFW_KEY_ENTER) {
            change_mode_active = false;
            printf("\n\tENTER PRESSED, CHANGE MODE DISABLED");
        }
        else {
            if(key == GLFW_KEY_BACKSPACE && change_mode_buffer.length() > 0) {
                change_mode_buffer.erase(change_mode_buffer.length() - 1);
            } else {

                //this is stupid but i dont want to sanitize it so idfc
                change_mode_buffer.push_back(key);
                printf("\n\tKEY ADDED TO BUFFER");
            }
        }

        return;
    }
    if(key == GLFW_KEY_ENTER) {
        printf("\n\tENTER PRESSED, CHANGE MODE ACTIVE");
        change_mode_active = true;
        return;
    }

    if(key == GLFW_KEY_1) {
        request_scene_type = SCENE_TYPE_TANK;
    } else if(key == GLFW_KEY_2) {
        request_scene_type = SCENE_TYPE_WIND_TUNNEL;
    } else if(key == GLFW_KEY_3) {
        request_scene_type = SCENE_TYPE_PAINT;
    }

    if(key == GLFW_KEY_P) {
        request_show_pressure = !request_show_pressure;
        show_pressure_setting = request_show_pressure;
        if(request_show_pressure == true) {
            request_show_temperature = false;
        } else if(show_temperature_setting == true) {
            request_show_temperature = true;
        }
    }
    if(key == GLFW_KEY_V) {
        request_show_velocity = !request_show_velocity;
    }
    if(key == GLFW_KEY_S) {
        request_show_smoke = !request_show_smoke;
    }
    if(key == GLFW_KEY_L) {
        request_show_streamlines = !request_show_streamlines;
    }
    if(key == GLFW_KEY_T) {
        request_show_temperature = !request_show_temperature;
        show_temperature_setting = request_show_temperature;
        if(request_show_temperature == true) {
            request_show_pressure = false;
        } else if(show_pressure_setting == true) {
            request_show_pressure = true;
        }
    }
    if(key == GLFW_KEY_C) {
        request_clear_pressure = !request_clear_pressure;
    }

    if(key == GLFW_KEY_MINUS) {
        request_delta_pressure_iterations -= 10;
    }
    if(key == GLFW_KEY_EQUAL) {
        request_delta_pressure_iterations += 10;
    }

    if(key == GLFW_KEY_SPACE) {
        paused = !paused;

    }
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0,0, width, height);
}

void terminate(GLFWwindow* window) {
    if(window != nullptr) {
        glfwDestroyWindow(window);
    }
    glfwTerminate();
}

void set_pixel_buffer(Scene* scene, unsigned char* buffer) {
    State* s = scene->s;

    fluid_type min_pressure = 1000000.0;
    fluid_type max_pressure = -10000.0;

    fluid_type min_temp = 100000000.0;
    fluid_type max_temp = -100000000.0;

    int nonzero_touched_values[4];
    int id_by_touched_value[4];

    int touched_values = 0;
    for(int i = 0; i < 4; i++) {
        nonzero_touched_values[i] = 0;
    }

    for(int ti = 0; ti < s->num_cells; ti++) {
        int touch_value = s->debug_touch_buffer[ti];
        if(touch_value == 0) {
            continue;
        }

        for(int existing_ti = 0; existing_ti < 4; existing_ti++) {
            if(nonzero_touched_values[existing_ti] == 0) {
                nonzero_touched_values[existing_ti] = touch_value;
                touched_values++;
                break;
            }
            if(nonzero_touched_values[existing_ti] == touch_value) {
                break;
            }
        }
    }

    unsigned char debug_colors[3 * 4] = {255,255,255, 50,50,50, 150,50,50, 50,150,50};

    for(int i = 0; i < s->num_cells; i++) {
        if(s->pressure[i] < min_pressure) {
            min_pressure = s->pressure[i];
        }
        if(s->pressure[i] > max_pressure) {
            max_pressure = s->pressure[i];
        }

        if(s->temperature_current[i] < min_temp) {
            min_temp = s->temperature_current[i];
        }
        if(s->temperature_current[i] > max_temp) {
            max_temp = s->temperature_current[i];
        }
    }

    int fluid_cells = 0;
    int solid_cells = 0;

    for(int y = 0; y < s->cells[1]; y++) {
        for(int x = 0; x < s->cells[0]; x++) {
            unsigned char color[4] = {255,255,255,255};
            int us = x + s->cells[0] * y;

            if(scene->show_pressure == true) {
                fluid_type pressure = s->pressure[us];
                getSciColor(color, pressure, min_pressure, max_pressure);

                if(scene->show_smoke == true) {
                    fluid_type sol = s->m_current[us];
                    /*if(sol == 1.0) {
                        color[0] = 0;
                        color[1] = 0;
                        color[2] = 0;
                    } *///else {

                    color[0] = (unsigned char)max(0.0, (fluid_type)color[0] - 255.0 * sol);
                    color[1] = (unsigned char)max(0.0, (fluid_type)color[1] - 255.0 * sol);
                    color[2] = (unsigned char)max(0.0, (fluid_type)color[2] - 255.0 * sol);
                    //}
                }
            }

            else if(scene->show_temperature == true) {
                fluid_type temp = s->temperature_current[us];
                getSciColor(color, temp, min_temp, max_temp);

                if(scene->show_smoke == true) {
                    fluid_type smoke = s->m_current[us];


                    color[0] = (unsigned char)max(0.0, (fluid_type)color[0] - 255.0 * smoke);
                    color[1] = (unsigned char)max(0.0, (fluid_type)color[1] - 255.0 * smoke);
                    color[2] = (unsigned char)max(0.0, (fluid_type)color[2] - 255.0 * smoke);
                }
            }

            else if(scene->show_smoke == true) {
                fluid_type sol = s->m_current[us];
                color[0] = (unsigned char)(255.0 * sol);
                color[1] = (unsigned char)(255.0 * sol);
                color[2] = (unsigned char)(255.0 * sol);
            }
            if(s->solids[x + s->cells[0] * y] == S_FIELD_FLUID) {
                //printf("\n\t\tFLUID AT %i, %i", x, y);
                fluid_cells++;
            }
            if(s->solids[x + s->cells[0] * y] == S_FIELD_SOLID) {
                solid_cells++;
                //if(x != 0 && x != s->cells[0] - 1 && y != 0 && y != s->cells[1] - 1) {

                //    printf("\n\t\tSOLID AT %i %i", x, y);
                //}
                color[0] = 30;
                color[1] = 30;
                color[2] = 30;
            }

            int red_x_start = 100;
            int red_y_start = 150;
            int red_1_width = 100;
            int red_1_height = 3;

            int red_x2_start = red_x_start + red_1_width - 1 - 3;
            int red_y2_start = red_y_start;
            int red_2_width = 3;
            int red_2_height = 50;

            //color[0] = (unsigned char) (x * 255 / s->cells[0]);
            //color[1] = (unsigned char) (y * 255 / s->cells[1]);
            //color[2] = 0;

            /*
            if(x >= red_x_start && x <= red_x_start + red_1_width - 1 && y >= red_y_start && y <= red_y_start + red_1_height - 1) {
                color[0] = 0;
                color[1] = 255;
                color[2] = 0;
            }
            if(x >= red_x2_start && x <= red_x2_start + red_2_width - 1 && y >= red_y2_start && y <= red_y2_start + red_2_height - 1) {
                color[0] = 0;
                color[1] = 0;
                color[2] = 255;
            }

            if(x == (int)s->cells[0] - 1){ //(s->cells[0] / 2)
                color[0] = 255;
                color[1] = 255;
                color[2] = 255;
            }
            */

#ifdef WHITEOUT_UNTOUCHED_CELLS
            if(s->debug_touch_buffer[x + s->cells[0] * y] != 0) {
                for(int ti = 0; ti < 4; ti++) {
                    if(nonzero_touched_values[ti] == s->debug_touch_buffer[x + s->cells[0] * y]) {
                        color[0] = debug_colors[3 * ti];
                        color[1] = debug_colors[3 * ti + 1];
                        color[2] = debug_colors[3 * ti + 2];
                        break;
                    }
                }
            }
#endif

            buffer[0 + 4 * (x + s->cells[0] * y)] = color[0];
            buffer[1 + 4 * (x + s->cells[0] * y)] = color[1];
            buffer[2 + 4 * (x + s->cells[0] * y)] = color[2];
            buffer[3 + 4 * (x + s->cells[0] * y)] = 255;
        }
    }


    //printf("\n----------------------------%i FLUID, %i SOLID------------------------------", fluid_cells, solid_cells);

}

unsigned char* resize_pixel_buffer(unsigned char* old_buffer, int old_width, int old_height, int new_width, int new_height) {
    unsigned char* new_buffer = new unsigned char[new_width * new_height * 4];
    delete [] old_buffer;
    for(int i = 0; i < new_width * new_height * 4; i++) {
        new_buffer[i] = 0;
    } 
    return new_buffer;
}

struct gForce {

};

struct ParsedNum {
    double num;
    unsigned int start_index;
    unsigned int end_index;
};


struct LogRelevantSimParams {
    public: 
        string scenario_ident;
        fluid_type dt;
        fluid_type h;
        int num_divergence_iterations;
        fluid_type density;
        int num_cells;
        fluid_type viscosity;
        fluid_type over_relax;

        static LogRelevantSimParams init(Scene* scene, fluid_type dt) {
            State* s = scene->s;
            string ident = "";
            switch(scene->scene_type) {
                case SCENE_TYPE_TANK:
                    ident = "Tank";
                    break;
                case SCENE_TYPE_WIND_TUNNEL:
                    ident = "Wind_Tunnel";
                    break;
                case SCENE_TYPE_PAINT:
                    ident = "Paint";
                    break;
            }
            return LogRelevantSimParams(
                ident, 
                dt,
                s->h,
                s->num_divergence_iterations,
                s->density,
                s->num_cells,
                s->viscosity,
                s->over_relaxation);
        }

        string to_string() {
            string ret = scenario_ident + "-";
            ret += std::to_string(dt) + "-";
            ret += std::to_string(h) + "-";
            ret += std::to_string(num_divergence_iterations) + "-";
            ret += std::to_string(density) + "-";
            ret += std::to_string(num_cells) + "-";
            ret += std::to_string(viscosity) + "-";
            ret += std::to_string(over_relax);
            return ret;
        }

        static optional<LogRelevantSimParams> from_string(string input) {
            stringstream stream(input);
            LogRelevantSimParams ret("", 0, 0, 0, 0, 0, 0, 0);
            string state = input;
            auto pos = state.find("-");
            int nth = 0;
            while(pos != string::npos && nth < 8) {
                string next = state.substr(0, pos);
                state.erase(0, pos + 1);

                switch(nth) {
                    case 0:
                        ret.scenario_ident = next;
                        break;
                    case 1:
                        ret.dt = (fluid_type)stod(next);
                        break;
                    case 2:
                        ret.h = (fluid_type)stod(next);
                        break;
                    case 3: 
                        ret.num_divergence_iterations = (int)stoi(next);
                        break;
                    case 4:
                    
                        ret.density = (fluid_type)stod(next);
                        break;
                    case 5:
                    
                        ret.num_cells = (int)stoi(next);
                        break;
                    case 6:
                    
                        ret.viscosity = (fluid_type)stod(next);
                        break;
                    case 7:
                        ret.over_relax = (fluid_type)stod(next);
                        break;
                }

                pos = state.find("-");
                nth++;
            }
            if(nth < 8) {
                return {};
            }
            return ret;
        }

        bool is_equal(LogRelevantSimParams* other) {
            return scenario_ident == other->scenario_ident
                && dt == other->dt
                && h == other->h
                && num_divergence_iterations == other->num_divergence_iterations
                && density == other->density
                && num_cells == other->num_cells
                && viscosity == other->viscosity
                && over_relax == other->over_relax;
        }

    private:
        LogRelevantSimParams(
            string scenario, 
            fluid_type dt,
            fluid_type h,
            int num_divergence_iterations,
            fluid_type density,
            int num_cells,
            fluid_type viscosity,
            fluid_type over_relax
            ) : 
            scenario_ident(scenario),
            dt(dt),
            h(h),
            num_divergence_iterations(num_divergence_iterations),
            density(density),
            num_cells(num_cells),
            viscosity(viscosity),
            over_relax(over_relax) 
            {}

    
};


string get_next_log_file_title(LogRelevantSimParams* params, filesystem::path logs_folder) {

    string title_stub = params->to_string();
    int num_matching = 0;

    for(const auto& log_file_entry : filesystem::recursive_directory_iterator(logs_folder)) {
        string log_file_name = log_file_entry.path().stem().string();
        auto returned = LogRelevantSimParams::from_string(log_file_name);
        if(!returned) {
            continue;
        }
        LogRelevantSimParams file_params = returned.value();
        if(params->is_equal(&file_params)) {
            num_matching++;
        }
    }

    return title_stub + "-" + to_string(num_matching) + ".csv";
}

int main() {
    if(!glfwInit()) {
        printf("\ninit failed!");
        return 0;
    }
    glfwSetErrorCallback(error_callback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    int window_width = INITIAL_WINDOW_WIDTH;
    int window_height = INITIAL_WINDOW_HEIGHT;

    int initial_window_width = window_width;
    int initial_window_height = window_height;

    GLFWwindow* window = glfwCreateWindow(window_width, window_height, "Fluid Simulator", NULL, NULL);
    if(!window) {
        printf("\nwindow initialization failed!");
        terminate(nullptr);
        return -1;
    }

    printf("\nwindow initialization succeeded!");

    glfwMakeContextCurrent(window);

    if(!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        printf("\nGLAD did not initialize!");
        terminate(window);
        return -1;
    }


    printf("\nglad initialization succeeded!");


    glViewport(0,0, window_width, window_height);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    glfwSetKeyCallback(window, key_callback);


    glfwSwapInterval(1);

    const char *vertex_shader_source = "#version 430 core\n"
        "layout (location = 0) in vec3 aPos;\n"
        "layout (location = 1) in vec2 tex_coord;\n"
        "uniform float time;\n"
        "out vec2 tex_coord_0;"
        "void main()\n"
        "{\n"
        "   tex_coord_0 = tex_coord;\n"
        "   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);\n"
        "}\0";


    unsigned int vertex_shader;
    vertex_shader = glCreateShader(GL_VERTEX_SHADER);

    glShaderSource(vertex_shader, 1, &vertex_shader_source, NULL);
    glCompileShader(vertex_shader);

    const char *fragment_shader_source = "#version 430 core\n"
        "in vec2 tex_coord_0;\n"
        "out vec4 FragColor;\n"
        "uniform sampler2D texture1;\n"
        "void main()\n"
        "{\n"
        "   FragColor = texture(texture1, tex_coord_0);\n"
        "}\n\0"; 

    ////vec4(1.0f, 0.5f, 0.2f, 1.0f);

    unsigned int fragment_shader;
    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    unsigned int shader_program;
    shader_program = glCreateProgram();


    const char *incompressability_compute_shader_source = R"(
        #version 430

        layout(local_size_x = 16, local_size_y = 16) in;

        layout(binding = 0) buffer Solids {
            float solids[];
        };

        layout(binding = 1) buffer UCurrent {
            float u_current[];
        };

        layout(binding = 2) buffer VCurrent {
            float v_current[];
        };

        layout(binding = 3) buffer Pressure {
            float pressure[];
        };

        uniform int n;
        uniform float cp;
        uniform float over_relaxation;
        uniform float viscosity;

        void main() {
            ivec2 gid = ivec2(gl_GlobalInvocationID.xy);
            int x = gid.x;
            int y = gid.y;
            int us = y * n + x;

            if (x < 1 || x >= n - 1 || y < 1 || y >= n - 1 || solids[us] == 1.0) {
                return;
            }

            int left = us - 1;
            int right = us + 1;
            int up = (y + 1) * n + x;
            int down = (y - 1) * n + x;

            float sol = solids[us];
            float sx0 = solids[left];
            float sx1 = solids[right];
            float sy0 = solids[down];
            float sy1 = solids[up];
            sol = sx0 + sx1 + sy0 + sy1;

            if (sol == 1.0) {
                return;
            }

            float divergence = u_current[right] - u_current[us] + v_current[up] - v_current[us];
            float p = -divergence / sol;
            p *= over_relaxation;
            pressure[us] += cp * p;

            u_current[us] -= sx0 * p / viscosity;
            u_current[right] += sx1 * p / viscosity;
            v_current[us] -= sy0 * p / viscosity;
            v_current[up] += sy1 * p / viscosity;
        }
        )";
    unsigned int compute_shader;
    compute_shader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(compute_shader, 1, &incompressability_compute_shader_source, NULL);
    glCompileShader(compute_shader);


    unsigned int compute_program;
    compute_program = glCreateProgram();
    glAttachShader(compute_program, compute_shader);
    glLinkProgram(compute_program);

    int success;
    char info_log[512];

    glGetProgramiv(compute_program, GL_LINK_STATUS, &success);
    if(success == 0) {
        glGetProgramInfoLog(compute_program, 512, NULL, info_log);
        cout << info_log << endl;
    }
    glUseProgram(compute_program);
    glDeleteShader(compute_shader);

    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    glLinkProgram(shader_program);

    int time_uniform_location = glGetUniformLocation(shader_program, "time"); 

    glGetProgramiv(shader_program, GL_LINK_STATUS, &success);
    if(success == 0) {
        glGetProgramInfoLog(shader_program, 512, NULL, info_log);
        cout << info_log << endl;
    }

    glUseProgram(shader_program);
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    //glDeleteShader(compute_shader);
    
    float vertices[] = {
            // positions          // colors           // texture coords
        1.0f,  -1.0f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // bottom right
        -1.0f, 1.0f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // top left
        -1.0f, -1.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
        1.0f,  1.0f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top right
    }; 
    

    /*
    float vertices[] = {
        // positions          // colors           // texture coords
         0.5f,  0.5f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // top right
         0.5f, -0.5f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // bottom right
        -0.5f, -0.5f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // bottom left
        -0.5f,  0.5f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // top left 
    };
    */

    unsigned int indices[] = {  
        0, 1, 2, // first triangle
        1, 0, 3  // second triangle
    };
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    // -------------------------
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
    // set the texture wrapping parameters
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);	// set texture wrapping to GL_REPEAT (default wrapping method)
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glEnable(GL_POLYGON_SMOOTH);



    //glGenerateMipmap(GL_TEXTURE_2D);

    Scene* scene = create_scene(STARTING_SCENARIO);

    State* s = scene->s;

    int color_channels = 4;
    unsigned char* cpu_pixel_buffer = new unsigned char[s->cells[1] * s->cells[0] * color_channels];
    for(int y = 0; y < s->cells[1]; y++) {
        for(int x = 0; x < s->cells[0]; x++) {
            int non_color_index = y * s->cells[0] * color_channels + x * color_channels;
            //R, G, B, A
            cpu_pixel_buffer[0 + non_color_index] = (unsigned char)255;
            cpu_pixel_buffer[1 + non_color_index] = (unsigned char)255;
            cpu_pixel_buffer[2 + non_color_index] = (unsigned char)255;
            cpu_pixel_buffer[3 + non_color_index] = (unsigned char)255;
        }
    }


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, s->cells[0], s->cells[1], 0, GL_RGBA, GL_UNSIGNED_BYTE, cpu_pixel_buffer);
    glGenerateMipmap(GL_TEXTURE_2D);

    int num_cells = s->num_cells;
    fluid_type dt = STARTING_DT;

    //compute shader stuff
    GLuint buffers[4];
    glGenBuffers(4,buffers);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[0]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, num_cells * sizeof(float), s->solids, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, buffers[0]);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[1]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, num_cells * sizeof(float), s->u_current, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, buffers[1]);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[2]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, num_cells * sizeof(float), s->v_current, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, buffers[2]);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[3]);
    glBufferData(GL_SHADER_STORAGE_BUFFER, num_cells * sizeof(float), s->pressure, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, buffers[3]);

    // Set uniform values
    glUniform1i(glGetUniformLocation(compute_program, "n"), s->cells[0]);
    glUniform1f(glGetUniformLocation(compute_program, "cp"), s->density * s->h / dt);
    glUniform1f(glGetUniformLocation(compute_program, "over_relaxation"), s->over_relaxation);
    glUniform1f(glGetUniformLocation(compute_program, "viscosity"), s->viscosity);

    // Dispatch the compute shader
    //glDispatchCompute((GLuint)n / 16, (GLuint)n / 16, 1);
    //glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    glUniform1i(glGetUniformLocation(shader_program, "texture1"), 0);

    int iterations = 0;
    fluid_type time = 0.0;
    fluid_type write_log_time = 10.0;
    fluid_type max_time = 10000.0;

    auto cur_path = filesystem::current_path();
    cout << endl << cur_path << endl;
    auto logs_folder = cur_path / "logs";

    cout << endl << logs_folder << endl;

    LogRelevantSimParams sim_params = LogRelevantSimParams::init(scene, dt);

    string log_file_title = get_next_log_file_title(&sim_params, logs_folder);
    ofstream log_file(logs_folder / log_file_title);
    log_file << "iteration,time,avg pressure,min pressure,max pressure,avg u velocity,avg v velocity,total energy\n";
    //[scenario ident]-[dt]-[h]-[num_divergence_iterations]-[density]-[num_cells]-[viscosity]-[over relax]-[number].csv
    // have a function that creates a stub file title given the scene and dt and then finds identical files in the logs folder
    // and sets the final log file title to number of matching stub names in logs folder + 1.

    //the file is not written if we dont reach an end state (total time elapsed >= write_log_time, or we close before then and when prompted say Y do write the log)
    //csv is of the form: iteration, time, avg pressure, min pressure, max pressure, avg u velocity, avg v velocity, total energy

    bool has_auto_paused_yet = false;
    bool is_auto_paused = false;

    printf("\n-----now starting main loop!-----");

    while(true) {
        clock_t loop_start = clock();
        if(is_sim_done(s, iterations, time, max_time)) {
            printf("\nmax time reached, returning");

            log_file.close();
            break;
        }

        if(!has_auto_paused_yet && PAUSE_AFTER > 0.0 && time >= PAUSE_AFTER) {
            paused = true;
            has_auto_paused_yet = true;
        }
        LogRow next_row;
        next_row.iteration = iterations;
        next_row.time = time;

        int old_window_width = window_width;
        int old_window_height = window_height;

        glfwGetFramebufferSize(window, &window_width, &window_height);
        glViewport(0,0, window_width, window_height);
        
        int old_scene_type = request_scene_type;
        bool old_change_mode_active = change_mode_active;

        for(int i = 0; i < 10; i++) {

            glfwPollEvents();
        }

        if(glfwWindowShouldClose(window)) {
            printf("\nwindow closed, returning");

            log_file.close();
            break;
        }

        if(change_mode_buffer.length() > 0) {

            //cout << change_mode_buffer << endl; 
            printf("\n\t\t CMB: %s", change_mode_buffer.c_str());
        }



        if(old_change_mode_active == true && change_mode_active == false) {
            printf("\nchange mode no longer active!");
            bool valid = false;
            bool found_command = false;
            for(int i = 0; i < commands.size(); i++) {
                if(change_mode_buffer.find(commands[i]) != string::npos) {
                    found_command = true;
                    break;
                }
            }
            if(found_command == true) {

                bool found_num = false;
                bool negative = false;
                bool less_than_one = false;
                int power = 0;
                double num = 0.0;
                //unsigned char[12] filter = "0123456789-.".c_str();
                const char* cmb = change_mode_buffer.c_str();
                vector<ParsedNum> nums;
                
                for(int i = 0; i < change_mode_buffer.length(); i++) {
                    if((cmb[i] < 48 || cmb[i] > 57) && cmb[i] != '-' && cmb[i] != '.') {
                        if(found_num == true) {
                            break;
                        } else {
                            if(negative == true) {
                                negative = false;
                            }
                            if(less_than_one == true) {
                                less_than_one = false;
                            }
                            continue;
                        }
                    }
                    if(found_num == false) {
                        if(cmb[i] == '-') {
                            if(negative == true) {
                                continue;
                            }
                            negative = true;
                            continue;
                        }
                        if(cmb[i] == '.') {
                            if(negative == true) {
                                negative = false;
                            }
                            less_than_one = true;
                            continue;
                        }
                        if(cmb[i] >= 48 && cmb[i] <= 57) { //i fucking love string processing
                            found_num = true;
                            valid = true;
                            ParsedNum n;
                            n.start_index = i;
                            nums.push_back(n);
                        }
                    }
                    if(found_num == true) {
                        if(cmb[i] == '-') {
                            nums[nums.size() - 1].end_index = i - 1;
                            nums[nums.size() - 1].num = num;
                            break;
                        }
                        if(cmb[i] == '.') {
                            less_than_one = true;
                            continue;
                        }
                        if(cmb[i] >= 48 && cmb[i] <= 57) {
                            int digit = cmb[i] - 48;
                            if(negative) {
                                digit *= -1;
                            }
                            if(less_than_one) {
                                num += pow(10, -power) * digit;
                                power++;
                            } else {
                                num *= 10.0;
                                num += (double)digit;
                            }
                        }
                    }
                }

                if(valid) {
                    if(change_mode_buffer.find(COMMAND_DIVERGENCE) != string::npos) {
                        s->num_divergence_iterations = (int)num;
                    } else if(change_mode_buffer.find(COMMAND_VISCOSITY) != string::npos) {
                        s->viscosity = num;
                    } else if(change_mode_buffer.find(COMMAND_PULSE) != string::npos) {
                        
                        for(int y = 0; y < s->cells[1]; y++) {
                            for(int x = 0; x < s->cells[0]; x++) {

                                if(x == 1) {
                                    s->u_current[x + s->cells[0] * y] = num;
                                    //s->v_current[x + s->cells[0] * y] = in_velocity / 0.50;
                                    //s->u_source[x + s->cells[0] * y] = in_velocity;
                                }
                                //if(x == (s->cells[0] - 3)) {
                                //    s->u_current[x + s->cells[0] * y] = 1.0 * in_velocity;
                                //}
                                //if(x == s->cells[0] - 1) {
                                //    s->u_current[x + s->cells[0] * y] = -in_velocity;
                                //}
                            }
                        }
                    } else if(change_mode_buffer.find(COMMAND_GRAV) != string::npos) {
                        s->g_strength = num;
                    } else if(change_mode_buffer.find(COMMAND_PRINT_PERF) != string::npos) {
                        print_perf = !print_perf;
                    }
                }
            }
            change_mode_buffer.erase();
        }

        if(request_scene_type != old_scene_type) {
            string msg;
            bool proceed = true;
            if(request_scene_type == SCENE_TYPE_PAINT) {
                msg = string("paint");
            }
            else if(request_scene_type == SCENE_TYPE_TANK) {
                msg = string("tank");
            }
            else if(request_scene_type == SCENE_TYPE_WIND_TUNNEL) {
                msg = string("wind tunnel");
            }
            else {
                proceed = false;
            }
            if(proceed == true) {
                printf("\n----------Scene type was changed to %s!----------", msg);
                scene->setup_scene(request_scene_type);


                log_file.close();
                if(time < MIN_TIME_TO_LOG) {
                    remove((logs_folder / log_file_title).u8string().c_str());
                }

                time = 0;
                iterations = 0;
                next_row.iteration = 0;
                next_row.time = 0;


                sim_params = LogRelevantSimParams::init(scene, dt);
                log_file_title = get_next_log_file_title(&sim_params, logs_folder);
                log_file = ofstream(logs_folder / log_file_title);

                log_file << "iteration, time, avg pressure, min pressure, max pressure, avg u velocity, avg v velocity, total energy\n";

                printf("\n\tNOW WRITING TO %s", log_file_title.c_str());
            }
        }
        scene->show_pressure = request_show_pressure;
        scene->show_smoke = request_show_smoke;
        scene->show_streamlines = request_show_streamlines;
        scene->show_velocity = request_show_velocity;
        scene->show_temperature = request_show_temperature;
        s->clear_pressure = request_clear_pressure;
        s->num_divergence_iterations += request_delta_pressure_iterations;
        request_delta_pressure_iterations = 0;
        fluid_type last_tick_time = 0.0;

        if(paused == false) {

            clock_t sim_start = clock();
            sim_step(s, iterations, time, dt, &next_row, scene->scene_type, 16, shader_program, compute_program);
            clock_t sim_end = clock();

            next_row.write(log_file);
            last_tick_time = (fluid_type)(sim_end - sim_start) / (CLOCKS_PER_SEC / 1000.0);
        }

        clock_t buffer_start = clock();
        set_pixel_buffer(scene, cpu_pixel_buffer);
        clock_t buffer_end = clock();

        fluid_type avg_divergence = calc_avg_divergence(s);
        if(!paused) {
            printf("\n\tavg divergence is %g", avg_divergence);
        }

        string scenario = "";
        if(scene->scene_type == SCENE_TYPE_TANK) {
            scenario = "Tank";
        }
        if(scene->scene_type == SCENE_TYPE_WIND_TUNNEL) {
            scenario = "Wind Tunnel";
        }
        if(scene->scene_type == SCENE_TYPE_PAINT) {
            scenario = "Paint";
        }
        string title = scenario + + ": iterations " + to_string(iterations) + ", time " + to_string(time) + "s, dt " + to_string(dt) + ", h: " + to_string(s->h) + " cells_x: " + to_string(s->cells[0]) + " cells_y: " + to_string(s->cells[1]);
        title += ", last tick time: " + to_string(last_tick_time) + ", over relax: " + to_string(s->over_relaxation) + ", divergence iters: " + to_string(s->num_divergence_iterations);
        title += ", avg divergence: " + to_string(avg_divergence);
        //format("{}: iteration {} dt {} h {} cells_x: {} cells_y: {}", scenario, iterations, dt, s->h, s->cells[0], s->cells[1]);
        glfwSetWindowTitle(window, title.c_str());
        //printf("\n\tset_pixel_buffer() took %g milliseconds", (fluid_type)(buffer_end - buffer_start) / (CLOCKS_PER_SEC / 1000.0));

        clock_t opengl_start = clock();


        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, s->cells[0], s->cells[1], 0, GL_RGBA, GL_UNSIGNED_BYTE, cpu_pixel_buffer);
        glGenerateMipmap(GL_TEXTURE_2D);


        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);	// set texture wrapping to GL_REPEAT (default wrapping method)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glBindTexture(GL_TEXTURE_2D, texture);
        glUseProgram(shader_program);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_TRIANGLES, 0, 6);

        glfwSwapBuffers(window);

        clock_t opengl_end = clock();

        //printf("\n\topengl graphics calls took %g milliseconds", (fluid_type)(opengl_end - opengl_start) / (CLOCKS_PER_SEC / 1000.0));

        clock_t loop_end = clock();


        if(paused == false) {
            printf("\niteration %i & time %g complete in %g milliseconds", iterations, time, (fluid_type)(loop_end - loop_start) / (CLOCKS_PER_SEC / 1000.0));

            time += dt;
            iterations++;
        }
    }

    terminate(window);
    return -1;
}

/*
int main2() {
    printf("\ncreating window");

    

    Window* window = new Window();

    bool running = true;
    while(running) {
        if(window->ProcessMessages() == false) {
            printf("\nrunning = false");
            running = false;
        }

        Sleep(10);
    }
    
    delete window;
    return 0;
}
*/
