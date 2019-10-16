#include "lj_compare.h"

real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path);
real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam,
	real sm_cunn);
real PRESSURE_DROP_SOURCE(real c_porosity, real vcsity,
	real ad_dynamic_shape_fac, real particle_density,
	real geometric_mean_diameter,
	real geometric_standard_deviation, real face_velocity,
	real face_dust_load);
real PRESSURE_DROP_SOURCE_DS(real c_porosity, real vcsity,
	real ad_dynamic_shape_fac, real particle_density,
	real geometric_mean_diameter,
	real geometric_standard_deviation, real face_velocity);
real TOTALL_porosity(real face_velocity, real face_dust_load);
real POROSITY_REAL(real dust_layer_height, real dust_loading_faceA, real dust_loading_mass, real dust_den);
real DUST_PERMEABLILITY(real dust_pressuredrop, real indust_velocity,
	real visco);
/*calculate for Re number based on pure FLUID*/
real RE_FLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_den, real characteristic_length);
/*calculate for Re number based on particle flow*/
real REN_PFLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_den, real par_vel[ND_ND], real par_diam);
/*Calculate for Re number based on particle flow with void ratio*/
real REP_PFLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_void, real fluid_den, real par_vel[ND_ND], real par_diam);
/*Theory derive from the paper "MP-PIC simulation of CFB riser with EMMS-based drag model" */
/*Egrun and wen yu*/
real COFF_B_P_MPPIC(real g_den, real g_void, real g_vel[ND_ND], real g_vis, real s_void, real p_vel[ND_ND], real p_diam, real coff_drag_p_mppic);
real COFF_DRAG_P_MPPIC(real REP_NUMBER);
/*smooth sphere particles in the normal flow. \
theory form the paper " An investigation of particle trajectories in two-phase flow systems"*/
real COFF_DRAG_P_SMOOTH(real REN_NUMBER);
/*for calculate the pressure gradient force in the MP-PIC *******/
real ACC_P_PG(real g_PG, real p_PG, real p_den, real p_void);
/*p_PG_star and Alfa are model parameters, p_void_max is s the solid volume fraction at close pack.*/
real ACC_P_PRESSURE(real Alfa, real p_PG_star, real p_void_max, real p_void);
/*gravity for the particle acceleration*/
real ACC_P_GVI(real g_den, real p_den, real g_acc, int dim);
/*distance in one time step*/
real POS_UPDATE(real pos0, real time_step, real p_acc, real p_vel0);
real DRAG_FLT(real drag_coff, real ren_number, real gas_vel, real par_vel, real gas_vis, real par_diam, real par_den);
real HIT_DIS_TIME(real pos0[ND_ND], real hit_pos[ND_ND], real vel0[ND_ND], real acc0[ND_ND], int dim);
/*theory derive from the paper "Effects of wall roughness on drag and lift forces of a particle at finite Reynolds number"
the rep_number should between 2 and 100.*/
real NEAR_DRAG_COFFI(real dis_towall, real rep_number);
/*theory derive from the paper "Effects of wall roughness on drag and lift forces of a particle at finite Reynolds number"
the rep_number should between 2 and 100.*/
real SAFF_LIFT_COFFI(real dis_towall, real rep_number);
int VORTI_VEC(real vorti_vec[ND_ND], cell_t c, Thread *t);
real REG_PFLUID(real vorti_mag, real gas_vis, real gas_den, real par_diam);
real CROSS_PF_VORTI(real cross_pf_vo[ND_ND], real par_v[ND_ND], real gas_v[ND_ND], real vort[ND_ND]);
real COMPOSITE_DRAG_COFFI(real dis_towall_coff, real rep_number);
int COMPOSITE_DRAG_FORCE(real drag_force[ND_ND], real gas_vel[ND_ND], real gas_den, real par_vel[ND_ND], real par_diam, real composite_drag_coff);
real COMPOSITE_LIFT_COFFI(real dis_towall, real rep_number, real par_diam);
real CORRELATION_LIFT_COFFI(real dis_towall_coff, real rep_number_G);
real DIS_TO_WALL_COFFI(real dis_towall, real par_diam);
real CORRELATION_LIFT_FORCE(real lift_force[ND_ND], real par_diam, real par_den, real par_vel[ND_ND], real gas_vis, real gas_den, real gas_vel[ND_ND], real vorti[ND_ND], real corr_lift_coffi);
real RADIAL_ANGLE(int tube_num, real cell_cord[ND_ND], real tube_center[ND_ND][TUBE_NUMBER]);
real DIS_TO_DUSTLAYER(real tube_co[ND_ND], real tube_r, real dustlayer_height, real particle_co[ND_ND], real mesh_below_surface);
real FLUENT_DRAG_RETURN(real dragF_count[ND_ND], real gas_vel[ND_ND], real par_vel[ND_ND], real gas_vis, real par_den, real par_diam);
real COMPOSTIE_DRAG_RETURN(real composite_drag_coff, real gas_den, real gas_vis, real par_den, real par_diam, real gas_vel[ND_ND], real par_vel[ND_ND]);
int LINE_TO_CYLINDER_INTERSECTION(real cylinder_o[ND_ND], real cy_height, real major_axis, real minor_axis, real pd_old[ND_ND], real pd_new[ND_ND], real line_cyc_hitpoint[ND_ND]);
int LINE_TO_CIRCLE_INTERSECTION(real cx, real cy, real r, real x1, real y1, real x2, real y2, real hitp[ND_ND]);
int ACC_TOTAL(real acc_tol[ND_ND], real dis_to_wall, real vorti[ND_ND], real gas_vis, real gas_v[ND_ND], real gas_den, real gas_p_g[ND_ND], real par_v[ND_ND], real par_den, real par_diam, real par_mass);
/*theory form the paper "Compression properties of dust cake of fine fly ashes
from a
fluidized bed coal combustion on a ceramic filter" and the formulation was chose
as e0. */
real POROUS_FROM_BOTTOM(real local_dust_height);
/*formulation for counting the dust-layer's height*/
real DUST_LOAD_HEIGHT(real dust_porosity, real dust_f_load, real particle_den);
int MEASSURE_TO_LAST(real* last_cord, real* meassured_cord, real distance_2, real prx_factor);
real LENGTH_TWO_POINTS(real* one, real* two, int dim);
/*using for outputting of dustlayer's thickness in one tube. the tube's circle point was preset. the dustlayer's thickness was output to a file including information
of tube-cell's location,cell's dustlayer height and cell's radial angel.
And in this code, the tube's axis parallel with Z direction. */
int DUSTLAYER_CELL_HEIGHT(int tube_num, Thread* t, cell_t c, real output_info[NEED_INFO]);
int NEAR_LAYER_NUM(cell_t c, Thread* t, int deposit_face_id);
/*count the face dust load in a cell's deposition face, mainly in x y direction
* , only for 3D*/
real FACE_DUST_LOAD_C(real total_mass_c, real deposition_face);