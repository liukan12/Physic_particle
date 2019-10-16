#include "lj_physics.h"

real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path)
{
	if (particle_diam > 2 * molecular_free_path) {
		return (1 + 2.468 * molecular_free_path / particle_diam);
	}
	else {
		return (1 + 3.294 * molecular_free_path / particle_diam);
	}
}

real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam,
	real sm_cunn)
{
	return (CUNNINGHAN_CORRECTION(vm_diam, MOLE_FREE_PATH) /
		CUNNINGHAN_CORRECTION(sm_diam, MOLE_FREE_PATH) * (vm_diam / sm_diam) *
		(vm_diam / sm_diam));
}

/*your pressure drop function*/
real PRESSURE_DROP_SOURCE(real c_porosity, real vcsity,
	real ad_dynamic_shape_fac, real particle_density,
	real geometric_mean_diameter,
	real geometric_standard_deviation, real face_velocity,
	real face_dust_load)
{
	real a_coefficient;
	real b_coefficient;
	real c_coefficient;
	real tpc_drop;
	c_coefficient =
		4 * log(geometric_standard_deviation) * log(geometric_standard_deviation);
	/*Message("your c_coefficient is %g\n", c_coefficient);*/
	a_coefficient = (2970 * vcsity * ad_dynamic_shape_fac * (1 - c_porosity)) /
		pow(c_porosity, 4);
	/*Message("your a_coefficient is %g\n", a_coefficient);*/
	b_coefficient = particle_density * geometric_mean_diameter * geometric_mean_diameter *
		exp(c_coefficient);
	/*Message("your b_coefficient is %g\n", b_coefficient);*/
	tpc_drop =
		a_coefficient * pow(b_coefficient, -1) * fabs(face_velocity) * fabs(face_dust_load);
	/*Message("your drop is %g\n", tpc_drop);*/
	return tpc_drop;
}

/*your pressure drop function*/
real PRESSURE_DROP_SOURCE_DS(real c_porosity, real vcsity,
	real ad_dynamic_shape_fac, real particle_density,
	real geometric_mean_diameter,
	real geometric_standard_deviation, real face_velocity)
{
	real a_coefficient;
	real b_coefficient;
	real c_coefficient;
	real tpc_drop;
	c_coefficient =
		4 * log(geometric_standard_deviation) * log(geometric_standard_deviation);
	/*Message("your c_coefficient is %g\n", c_coefficient);*/
	a_coefficient = (2970 * vcsity * ad_dynamic_shape_fac * (1 - c_porosity)) /
		pow(c_porosity, 4);
	/*Message("your a_coefficient is %g\n", a_coefficient);*/
	b_coefficient = particle_density * geometric_mean_diameter * geometric_mean_diameter *
		exp(c_coefficient);
	/*Message("your b_coefficient is %g\n", b_coefficient);*/
	tpc_drop =
		a_coefficient * pow(b_coefficient, -1)* fabs(face_velocity);
	/*Message("your drop is %g\n", tpc_drop);*/
	return tpc_drop;
}

real TOTALL_porosity(real face_velocity, real face_dust_load)
{
	return (1 -
		(0.88 -
			0.37 * pow(fabs(face_velocity), 0.36) *
			pow(fabs(face_dust_load), 0.25)) *
		pow(fabs(face_velocity), 0.36) * pow(fabs(face_dust_load), 0.25));
}

real POROSITY_REAL(real dust_layer_height, real dust_loading_faceA, real dust_loading_mass, real dust_den)
{
	return (dust_loading_mass / dust_den) / (dust_layer_height*dust_loading_faceA);
}

real DUST_PERMEABLILITY(real dust_pressuredrop, real indust_velocity,
	real visco)
{
	return (dust_pressuredrop / (indust_velocity * visco) * (-1));
}

/*calculate for Re number based on pure FLUID*/
real RE_FLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_den, real characteristic_length)
{
	real fluid_Vmag = 0.0;
#if RP_3D
	fluid_Vmag = sqrt(fluid_vel[0] * fluid_vel[0] + fluid_vel[1] * fluid_vel[1] + fluid_vel[2] * fluid_vel[2]);
#endif
	return (fluid_Vmag*fluid_den*characteristic_length / fluid_vis);
}

/*calculate for Re number based on particle flow*/
real REN_PFLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_den, real par_vel[ND_ND], real par_diam)
{
	return (DIS_3D(fluid_vel, par_vel)*fluid_den*par_diam / fluid_vis);
}

/*Calculate for Re number based on particle flow with void ratio*/
real REP_PFLUID(real fluid_vis, real fluid_vel[ND_ND], real fluid_void, real fluid_den, real par_vel[ND_ND], real par_diam)
{
	return (DIS_3D(fluid_vel, par_vel)*fluid_den*par_diam*fluid_void / fluid_vis);
}

/*Theory derive from the paper "MP-PIC simulation of CFB riser with EMMS-based drag model" */
/*Egrun and wen yu*/
real COFF_B_P_MPPIC(real g_den, real g_void, real g_vel[ND_ND], real g_vis, real s_void, real p_vel[ND_ND], real p_diam, real coff_drag_p_mppic)
{
	if (g_void < 0.8)
	{
		return (150 * g_vis*s_void / (g_void*p_diam*p_diam) + 1.75*g_den*DIS_3D(g_vel, p_vel) / p_diam);
	}
	else
	{
		return (0.75*g_den*g_void*DIS_3D(g_vel, p_vel) / p_diam * pow(g_void, -2.7)*coff_drag_p_mppic);
	}
}

real COFF_DRAG_P_MPPIC(real REP_NUMBER)
{
	if (REP_NUMBER < 1000.0)
	{
		return (24 * (1 + 0.15*pow(REP_NUMBER, 0.687)) / REP_NUMBER);
	}
	else
	{
		return (0.44);
	}
}

/*smooth sphere particles in the normal flow. \
theory form the paper " An investigation of particle trajectories in two-phase flow systems"*/
real COFF_DRAG_P_SMOOTH(real REN_NUMBER)
{
	real drag_coff = 0.;
	if (REN_NUMBER < 0.1)
	{
		drag_coff = 24.0 / REN_NUMBER;
	}
	else if (REN_NUMBER > 0.1&& REN_NUMBER < 1.0)
	{
		drag_coff = 22.73 / REN_NUMBER + 0.0903 / (REN_NUMBER*REN_NUMBER) + 3.69;
	}
	else if (REN_NUMBER > 1.0&& REN_NUMBER < 10.0)
	{
		drag_coff = 29.1667 / REN_NUMBER - 3.8889 / (REN_NUMBER*REN_NUMBER) + 1.222;
	}
	else if (REN_NUMBER > 10.0&& REN_NUMBER < 100.0)
	{
		drag_coff = 46.5 / REN_NUMBER - 116.67 / (REN_NUMBER*REN_NUMBER) + 0.6167;
	}
	else if (REN_NUMBER > 100.0&& REN_NUMBER < 1000.0)
	{
		drag_coff = 98.33 / REN_NUMBER - 2778.0 / (REN_NUMBER*REN_NUMBER) + 0.3644;
	}
	else if (REN_NUMBER > 1000.0&& REN_NUMBER < 5000.0)
	{
		drag_coff = 148.62 / REN_NUMBER - 47500.0 / (REN_NUMBER*REN_NUMBER) + 0.357;
	}
	else if (REN_NUMBER > 5000.0&& REN_NUMBER < 10000.0)
	{
		drag_coff = -490.546 / REN_NUMBER + 578700.0 / (REN_NUMBER*REN_NUMBER) + 0.46;
	}
	else if (REN_NUMBER > 10000.0&& REN_NUMBER < 50000.0)
	{
		drag_coff = -1662.5 / REN_NUMBER + 5416700.0 / (REN_NUMBER*REN_NUMBER) + 0.5191;
	}

	return drag_coff;
}

/*for calculate the pressure gradient force in the MP-PIC *******/
real ACC_P_PG(real g_PG, real p_PG, real p_den, real p_void)
{
	return (-(g_PG / p_den + p_PG / (p_void*p_den)));
}

/*p_PG_star and Alfa are model parameters, p_void_max is s the solid volume fraction at close pack.*/
real ACC_P_PRESSURE(real Alfa, real p_PG_star, real p_void_max, real p_void)
{
	return (p_PG_star*pow(p_void, Alfa) / (p_void_max - p_void));
}

/*gravity for the particle acceleration*/
real ACC_P_GVI(real g_den, real p_den, real g_acc, int dim)
{
	if (dim == G_DIRECTION)
		return (dim / abs(dim)*(p_den - g_den) / p_den * g_acc);
	else
		return 0;
}

/*distance in one step*/
real POS_UPDATE(real pos0, real time_step, real p_acc, real p_vel0)
{
	return (pos0 + time_step * p_vel0 + 0.5*time_step*time_step*p_acc);
}

real DRAG_FLT(real drag_coff, real ren_number, real gas_vel, real par_vel, real gas_vis, real par_diam, real par_den)
{
	return gas_vis / (par_den*par_diam*par_diam) * 18 * drag_coff*ren_number / 24 * (gas_vel - par_vel);
}

/*gravity,pressure gradient force,drag force and lift force*/
int ACC_TOTAL(real acc_tol[ND_ND], real dis_to_wall, real vorti[ND_ND], real gas_vis, real gas_v[ND_ND], real gas_den, real gas_p_g[ND_ND], real par_v[ND_ND], real par_den, real par_diam, real par_mass)
{
	real drag_coff = 0.0;
	int ii = 0, dim = 0;
	real coff_drag_p_mppic = 0.;
	real acc_drag[ND_ND] = { 0., };
	real acc_grav[ND_ND] = { 0., };
	real acc_lift[ND_ND] = { 0., };
	real acc_pres[ND_ND] = { 0., };
	real rep_number = 0.0;
	real dis_to_wall_coff = 0.;
	real comp_lift_coff = 0., comp_drag_coff = 0.;
	real lift_f[ND_ND], drag_f[ND_ND];

	rep_number = REN_PFLUID(gas_vis, gas_v, gas_den, par_v, par_diam);
	comp_lift_coff = COMPOSITE_LIFT_COFFI(dis_to_wall, rep_number, par_diam);
	CORRELATION_LIFT_FORCE(lift_f, par_diam, par_den, par_v, gas_vis, gas_den, gas_v, vorti, comp_lift_coff);
	dis_to_wall_coff = DIS_TO_WALL_COFFI(dis_to_wall, par_diam);
	/*coff_drag_p_mppic = COFF_DRAG_P_MPPIC(rep_number);*/
	/*coff_b_p = COFF_B_P_MPPIC(C_R(c, t), (1 - solid_prosity), fluid_v, AIR_VISOCOSITY, solid_prosity, par_v, TP_DIAM(tp), coff_drag_p_mppic);*/
	comp_drag_coff = COMPOSITE_DRAG_COFFI(dis_to_wall_coff, rep_number);
	COMPOSITE_DRAG_FORCE(drag_f, gas_v, gas_den, par_v, par_diam, comp_drag_coff);

	for (dim = 0; dim < ND_ND; dim++)
	{
		if (dim == abs(G_DIRECTION))
		{
			acc_grav[dim] = G_DIRECTION / abs(G_DIRECTION)*G_CONST*(par_den - gas_den) / par_den;
		}
		acc_pres[dim] = -gas_p_g[dim] / par_den;
		acc_drag[dim] = drag_f[dim];
		acc_lift[dim] = lift_f[dim] / par_mass;
		acc_tol[dim] = acc_grav[dim] + acc_pres[dim] + acc_drag[dim] + acc_lift[dim];
	}

	return 0;
}

real HIT_DIS_TIME(real pos0[ND_ND], real hit_pos[ND_ND], real vel0[ND_ND], real acc0[ND_ND], int dim)
{
	real dis_to_hit = 0., vel_mag = 0., acc_mag = 0., hit_time = 0.;
	if (dim == 2)
	{
		dis_to_hit = DIS_2D(pos0, hit_pos);
		vel_mag = MAG_2D(vel0);
		acc_mag = MAG_2D(acc0);
	}
	if (dim == 3)
	{
		dis_to_hit = DIS_3D(pos0, hit_pos);
		vel_mag = MAG_3D(vel0);
		acc_mag = MAG_3D(acc0);
	}
	if (acc_mag != 0.)
		hit_time = (sqrt(vel_mag*vel_mag + 2 * acc_mag*dis_to_hit) - vel_mag) / acc_mag;
	else
		hit_time = dis_to_hit / vel_mag;

	return hit_time;
}

/*theory derive from the paper "Effects of wall roughness on drag and lift forces of a particle at finite Reynolds number"
the rep_number should between 2 and 100.*/
real NEAR_DRAG_COFFI(real dis_towall, real rep_number)
{
	real Cds;
	real alfa_D;
	real beta_D;
	alfa_D = 0.15 - 0.046*(1.0 - 0.16*dis_towall*dis_towall)*exp(-0.7*dis_towall);
	beta_D = 0.687 + 0.066*(1.0 - 0.76*dis_towall*dis_towall)*exp(-pow(dis_towall, 0.9));
	Cds = 24.0 / rep_number * (1.0 + 0.138*exp(-2.0*dis_towall) + 9.0 / 16.0 / (1.0 + 2 * dis_towall))*(1.0 + alfa_D * pow(rep_number, beta_D));
	return Cds;
}

/*theory derive from the paper "Effects of wall roughness on drag and lift forces of a particle at finite Reynolds number"
the rep_number should between 2 and 100.*/
real SAFF_LIFT_COFFI(real dis_towall, real rep_number)
{
	real Cls;
	real alfa_L;
	real beta_L;
	real lamb_D;
	alfa_L = -exp(-0.3 + 0.025*rep_number);
	beta_L = 0.8 + 0.01*rep_number;
	lamb_D = (1.0 - exp(-dis_towall))*pow((rep_number / 250.0), 2.5);
	Cls = 3.663 / pow((rep_number*rep_number + 0.1173), 0.22)*exp(-0.5*dis_towall*pow((rep_number / 250.0), 4.0 / 3.0))*(exp(alfa_L*pow(dis_towall, beta_L)) - lamb_D);

	return Cls;
}

int VORTI_VEC(real vorti_vec[ND_ND], cell_t c, Thread *t)
{
	vorti_vec[0] = C_DWDY(c, t) - C_DVDZ(c, t);
	vorti_vec[1] = C_DUDZ(c, t) - C_DWDX(c, t);
	vorti_vec[2] = C_DVDX(c, t) - C_DUDY(c, t);

	return 0;
}

real REG_PFLUID(real vorti_mag, real gas_vis, real gas_den, real par_diam)
{
	return (gas_den*par_diam*par_diam*vorti_mag / gas_vis);
}

/* cross product of particle relative velocity and vorticity*/
real CROSS_PF_VORTI(real cross_pf_vo[ND_ND], real par_v[ND_ND], real gas_v[ND_ND], real vort[ND_ND])
{
	int dim;
	real rel_par_v[ND_ND];
	for (dim = 0; dim < ND_ND; dim++)
	{
		rel_par_v[dim] = gas_v[dim] - par_v[dim];
	}

	NV_CROSS(cross_pf_vo, rel_par_v, vort);

	return 0;
}

/*the theory derive from the paper
"The effects of near wall corrections to hydrodynamic forces on particle deposition and transport in vertical turbulent boundary layers"
on page 67 formula 23-24*/
real COMPOSITE_DRAG_COFFI(real dis_towall_coff, real rep_number)
{
	real Cd;
	real Cd0;

	Cd0 = (1.028 - 0.07 / (1 + 4 * dis_towall_coff*dis_towall_coff) - 8 / 15 * log((270 * dis_towall_coff / (135 + 256 * dis_towall_coff)))) * 24 / rep_number;
	Cd = (1.0 + 0.15*(1 - exp(-pow(dis_towall_coff, 0.5)))*pow(rep_number, (0.687 + 0.313*exp(-2 * pow(dis_towall_coff, 0.5)))))*Cd0;

	return Cd;
}

int COMPOSITE_DRAG_FORCE(real drag_force[ND_ND], real gas_vel[ND_ND], real gas_den, real par_vel[ND_ND], real par_diam, real composite_drag_coff)
{
	int ii;
	real contact_area = 0;
	contact_area = M_PI * par_diam*par_diam / 4;

	for (ii = 0; ii < ND_ND; ii++)
	{
		drag_force[ii] = -0.5*gas_den*DIS_3D(gas_vel, par_vel)*(gas_vel[ii] - par_vel[ii])*contact_area*composite_drag_coff;
	}

	return 0;
}

/*the theory derive from the paper
"The effects of near wall corrections to hydrodynamic forces on particle deposition and transport in vertical turbulent boundary layers"
on page 67 formula 33-39*/
real COMPOSITE_LIFT_COFFI(real dis_towall, real rep_number, real par_diam)
{
	real Clt;
	real Cltw;
	real Clt0 = 0.;
	real L0;
	real flr;
	real fr0;
	real fr1;
	real gr;
	real dis_towall_coff;

	dis_towall_coff = dis_towall / par_diam - 0.5;

	L0 = dis_towall * rep_number / par_diam;
	gr = 3 * exp(-0.17*pow(rep_number, 0.7));
	Cltw = 0.313 + 0.812*exp(-0.125*pow(rep_number, 0.77));
	if (L0 >= 0 && L0 < 10)
	{
		Clt0 = (9 / 8 + 5.78e-6*L0)*exp(-0.292*L0);
	}
	else if (L0 > 10)
	{
		Clt0 = 8.94*pow(L0, -2.09);
	}

	fr0 = 1 + 0.329*rep_number + 0.00485*rep_number*rep_number;
	fr1 = -0.9*tanh(0.022*rep_number);
	flr = fr0 * Clt0*pow(dis_towall, fr1);
	Clt = flr + (Cltw - fr0 * Clt0*pow(0.5, fr1))*exp(-11 * pow((dis_towall_coff / gr), 1.2));

	return Clt;
}

/*the theory derive from the paper
"The effects of near wall corrections to hydrodynamic forces on particle deposition and transport in vertical turbulent boundary layers"
on page 68 formula 31-32*/
real CORRELATION_LIFT_COFFI(real dis_towall_coff, real rep_number_G)
{
	real Cls;
	real Cls_W;
	real alfa_sL;
	real beta_sL;
	real lamb_sL;
	alfa_sL = -exp(-0.3 + 0.025*rep_number_G);
	beta_sL = 0.8 + 0.01*rep_number_G;
	lamb_sL = (1.0 - exp(-dis_towall_coff))*pow((rep_number_G / 250.0), 2.5);
	Cls_W = 3.663 / pow((rep_number_G*rep_number_G + 0.1173), 0.22);
	Cls = Cls_W * exp(-0.5*dis_towall_coff*pow((rep_number_G / 250.0), 4.0 / 3.0))*(exp(alfa_sL*pow(dis_towall_coff, beta_sL)) - lamb_sL);
	return Cls;
}

real DIS_TO_WALL_COFFI(real dis_towall, real par_diam)
{
	real dis_to_wall_coffi = dis_towall / par_diam - 0.5;
	if (dis_to_wall_coffi < MIN_COFF)
		dis_to_wall_coffi = MIN_COFF;

	return dis_to_wall_coffi;
}

real CORRELATION_LIFT_FORCE(real lift_force[ND_ND], real par_diam, real par_den, real par_vel[ND_ND], real gas_vis, real gas_den, real gas_vel[ND_ND], real vorti[ND_ND], real corr_lift_coffi)
{
	int ii;
	real cross_p_f[ND_ND] = { 0., }, contact_area = 0.;
	CROSS_PF_VORTI(cross_p_f, par_vel, gas_vel, vorti);
	contact_area = M_PI * par_diam*par_diam / 4;

	for (ii = 0; ii < ND_ND; ii++)
	{
		lift_force[ii] = 0.5*gas_den*contact_area*par_diam* cross_p_f[ii] * corr_lift_coffi;
	}

	return 0;
}
/*............................*/

real RADIAL_ANGLE(int tube_num, real cell_cord[ND_ND], real tube_center[ND_ND][TUBE_NUMBER])
{
	real radial_angle = 0.;
	real tube_co[ND_ND];
	int ii = 0;
	for (ii = 0; ii < ND_ND; ii++)
	{
		tube_co[ii] = tube_center[ii][tube_num];
	}

	if ((cell_cord[0] - tube_co[0]) >= 0)
	{
		if ((tube_co[0] == cell_cord[0]) && (cell_cord[1] - tube_co[1]) > 0)
			radial_angle = M_PI / 2;
		else if ((tube_co[0] == cell_cord[0]) && (cell_cord[1] - tube_co[1]) < 0)
			radial_angle = M_PI / 2 * 3;
		else
		{
			radial_angle = atan((cell_cord[1] - tube_co[1]) / (cell_cord[0] - tube_co[0]));
			if (radial_angle >= 0)
				radial_angle = radial_angle;
			else
				radial_angle = M_PI * 2 + radial_angle;
		}
	}
	else if ((tube_co[0] != cell_cord[0]) && (cell_cord[0] - tube_co[0]) < 0)
		radial_angle = M_PI + atan((cell_cord[1] - tube_co[1]) / (cell_cord[0] - tube_co[0]));

	return radial_angle;
}

real FLUENT_DRAG_RETURN(real dragF_count[ND_ND], real gas_vel[ND_ND], real par_vel[ND_ND], real gas_vis, real par_den, real par_diam)
{
	real dragcoff_fluent = 0.;
	dragcoff_fluent = fabs(dragF_count[0] / (gas_vel[0] - par_vel[0]) / (gas_vis / (par_den*par_diam*par_diam)));
	return dragcoff_fluent;
}

real COMPOSTIE_DRAG_RETURN(real composite_drag_coff, real gas_den, real gas_vis, real par_den, real par_diam, real gas_vel[ND_ND], real par_vel[ND_ND])
{
	real vel_mag = 0.;
	real contact_area = 0.;
	real fluent_drag_return = 0.;
	contact_area = par_diam * par_diam*M_PI / 4;
	vel_mag = DIS_3D(gas_vel, par_vel);
	fluent_drag_return = 0.5*gas_den*vel_mag*par_diam * par_diam*M_PI / 4 * composite_drag_coff / (gas_vis / (par_den*par_diam*par_diam));
	return fluent_drag_return;
}

/*the tube's axis parallel with Z direction*/
real DIS_TO_DUSTLAYER(real tube_co[ND_ND], real tube_r, real dustlayer_height, real particle_co[ND_ND], real mesh_below_surface)
{
	real pdis_to_tube = 0.;
	real dust_sur_to_tube = 0.;
	real dis_p_to_dust = 0.;
	pdis_to_tube = DIS_2D(particle_co, tube_co);
	dust_sur_to_tube = dustlayer_height + tube_r - mesh_below_surface;
	dis_p_to_dust = pdis_to_tube - dust_sur_to_tube;
	if (dis_p_to_dust < MIN_GAP)
		dis_p_to_dust = MIN_GAP;
	return dis_p_to_dust;
}

int LINE_TO_CYLINDER_INTERSECTION(real cylinder_o[ND_ND], real cy_height, real major_axis, real minor_axis, real pd_old[ND_ND], real pd_new[ND_ND], real line_cyc_hitpoint[ND_ND])
{
	int ii;
	real direction_vector[ND_ND] = { 0.,0.,0. };
	real coff_a = major_axis, coff_b = minor_axis; /*for the ellipse*/
	real x1[ND_ND] = { ERROR_POINTS,ERROR_POINTS ,ERROR_POINTS }, x2[ND_ND] = { ERROR_POINTS,ERROR_POINTS ,ERROR_POINTS };
	real dis_to_co = DIS_2D(pd_new, cylinder_o);
	real dis_1 = 0., dis_2 = 0.;
	real pd_old_lc[ND_ND] = { 0.,0.,0. };

	for (ii = 0; ii < ND_ND; ii++)
	{
		line_cyc_hitpoint[ii] = ERROR_POINTS;/*initialize the hitting point*/
		pd_old_lc[ii] = pd_old[ii] - cylinder_o[ii];/*convert into local coordination*/ 
	}

	if (dis_to_co > major_axis) /*if the particle have not moved into the cylinder zone*/
	{
		return 0;
	}
	else
	{
		/*line direction vector modified by x-axis*/

		if (pd_new[0] != pd_old[0])
		{
			direction_vector[0] = 1;
			direction_vector[1] = (pd_new[1] - pd_old[1]) / (pd_new[0] - pd_old[0]);
			direction_vector[2] = (pd_new[2] - pd_old[2]) / (pd_new[0] - pd_old[0]);

			/*delta>0--keep the real number */
			if ((pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old[0], 2) + 2 * direction_vector[1] * pd_old[0] * pd_old[1] - pow(pd_old[1], 2)) >= 0)
			{
				x1[0] = pd_old_lc[0] + (-pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] - pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2));
				x1[1] = direction_vector[1] * (-pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] - pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2)) + pd_old_lc[1];
				x1[2] = direction_vector[2] * (-pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] - pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2)) + pd_old_lc[2];
				x2[0] = pd_old_lc[0] - (pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] + pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2));
				x2[1] = -direction_vector[1] * (pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] + pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2)) + pd_old_lc[1];
				x2[2] = -direction_vector[2] * (pow(coff_a, 2)*direction_vector[1] * pd_old_lc[1] + pow(coff_b, 2)*pd_old_lc[0] + sqrt(pow(coff_a, 2)*pow(coff_b, 2)*(pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2) - pow(direction_vector[1], 2)*pow(pd_old_lc[0], 2) + 2 * direction_vector[1] * pd_old_lc[0] * pd_old_lc[1] - pow(pd_old_lc[1], 2)))) / (pow(coff_a, 2)*pow(direction_vector[1], 2) + pow(coff_b, 2)) + pd_old_lc[2];
				dis_1 = DIS_2D(pd_old_lc, x1);
				dis_2 = DIS_2D(pd_old_lc, x2);

				if (dis_1 <= dis_2) /*find nearest intersection and convert into global coordination*/
				{
					if (x1[2] <= cy_height)
					{
						line_cyc_hitpoint[0] = x1[0] + cylinder_o[0];
						line_cyc_hitpoint[1] = x1[1] + cylinder_o[1];
						line_cyc_hitpoint[2] = x1[2] + cylinder_o[2];
						return 1;
					}
					else
						return 0;
				}
				else
				{
					if (x2[2] <= cy_height)
					{
						line_cyc_hitpoint[0] = x2[0] + cylinder_o[0];
						line_cyc_hitpoint[1] = x2[1] + cylinder_o[1];
						line_cyc_hitpoint[2] = x2[2] + cylinder_o[2];
						return 1;
					}
					else
						return 0;
				}
			}
			else
				return 0; /*delta<0--drop the complex number*/
		}

		if (pd_new[0] == pd_old[0] && pd_new[1] != pd_old[1] && pd_new[2] != pd_old[2])
		{
			direction_vector[1] = 1;
			direction_vector[2] = (pd_new[2] - pd_old[2]) / (pd_new[1] - pd_old[1]);

			if (coff_a*coff_a - pd_old_lc[0] * pd_old_lc[0] >= 0)
			{
				x1[0] = pd_old_lc[0];
				x1[1] = sqrt((coff_a*coff_a - pd_old_lc[0] * pd_old_lc[0])*coff_b*coff_b / coff_a / coff_a);
				x1[2] = (x1[1] - pd_old_lc[1])*direction_vector[2] + pd_old_lc[2];
				x2[0] = pd_old_lc[0];
				x2[1] = -sqrt((coff_a*coff_a - pd_old_lc[0] * pd_old_lc[0])*coff_b*coff_b / coff_a / coff_a);
				x2[2] = (x2[1] - pd_old_lc[1])*direction_vector[2] + pd_old_lc[2];
				dis_1 = DIS_2D(pd_old_lc, x1);
				dis_2 = DIS_2D(pd_old_lc, x2);
				if (dis_1 <= dis_2) /*find nearest intersection and convert into global coordination*/
				{
					if (x1[2] <= cy_height)
					{
						line_cyc_hitpoint[0] = pd_old[0];
						line_cyc_hitpoint[1] = x1[1] + cylinder_o[1];
						line_cyc_hitpoint[2] = x1[2] + cylinder_o[2];
						return 1;
					}
					else
						return 0;
				}
				else
				{
					if (x2[2] <= cy_height)
					{
						line_cyc_hitpoint[0] = pd_old[0];
						line_cyc_hitpoint[1] = x2[1] + cylinder_o[1];
						line_cyc_hitpoint[2] = x2[2] + cylinder_o[2];
						return 1;
					}
					else
						return 0;
				}
			}
			else
				return 0;
		}

		if ((pd_new[0] == pd_old[0]) && (pd_new[1] == pd_old[1]))
		{
			if (pd_new[1] * pd_new[1] / (coff_b) / (coff_b)+pd_new[0] * pd_new[0] / (coff_a) / (coff_a) == 1)
			{
				line_cyc_hitpoint[0] = pd_new[0];
				line_cyc_hitpoint[1] = pd_new[1];
				line_cyc_hitpoint[2] = pd_new[2];
				return 1;
			}
			else
				return 0;
		}
		else
			return 0;
	}

	return 0;
}

int LINE_TO_CIRCLE_INTERSECTION(real cx, real cy, real r, real x1, real y1, real x2, real y2, real hitp[ND_ND])
{
	real dx = 0., dy = 0., dr = 0., D = 0., hitx1 = 0., hity1 = 0., hitx2 = 0., hity2 = 0.;
	real dis1 = 0., dis2 = 0., dis_to_co = 0.;
	real check_formula = 0.;
	real pos1[2] = { 0., }, center_o[2] = { 0., };
	pos1[0] = x2;
	pos1[1] = y2;
	center_o[0] = cx;
	center_o[1] = cy;
	dis_to_co = DIS_2D(pos1, center_o);

	if (dis_to_co > r)
	{
		return 0;
	}
	else
	{
		dx = x2 - x1;
		dy = y2 - y1;
		dr = sqrt(dx*dx + dy * dy);
		D = (x1 - cx)*(y2 - cy) - (x2 - cx)*(y1 - cy);

		check_formula = r * r*dr*dr - D * D;
		if (check_formula >= 0)
		{
			hitx1 = (D*dy + (sign(dy))*dx*sqrt(r*r*dr*dr - D * D)) / (dr*dr) + cx;
			hity1 = (-D * dx + fabs(dy)*sqrt(r*r*dr*dr - D * D)) / (dr*dr) + cy;
			hitx2 = (D*dy - (sign(dy))*dx*sqrt(r*r*dr*dr - D * D)) / (dr*dr) + cx;
			hity2 = (-D * dx - fabs(dy)*sqrt(r*r*dr*dr - D * D)) / (dr*dr) + cy;
			dis1 = (x1 - hitx1)*(x1 - hitx1) + (y1 - hity1)*(y1 - hity1);
			dis2 = (x1 - hitx2)*(x1 - hitx2) + (y1 - hity2)*(y1 - hity2);
			if (dis1 <= dis2)
			{
				hitp[0] = hitx1;
				hitp[1] = hity1;
			}
			else
			{
				hitp[0] = hitx2;
				hitp[1] = hity2;
			}
			return 1;
		}
		else
		{
			hitp[0] = ERROR_POINTS;
			hitp[1] = ERROR_POINTS;
			return 0;
		}

		return 0;
	}
}
/*theory form the paper "Compression properties of dust cake of fine fly ashes
from a
fluidized bed coal combustion on a ceramic filter" and the formulation was chose
as e0.*/
real POROUS_FROM_BOTTOM(real local_dust_height)
{
	real factor_a;
	real factor_b;
	real factor_c;
	real porous_in_cell = 0.92;
	factor_a = 0.81;
	factor_b = 0.13;
	factor_c = -8.374;
	porous_in_cell = factor_a +
		factor_b * exp(factor_c * (0.6 - (local_dust_height) * 1000));

	if (porous_in_cell > 1.0 || porous_in_cell < 0.0)
	{
		Message0("you get a wrong porous %g.\n", porous_in_cell);
		return 1.0;
	}
	else
	{
		/*if (DEBUG_TRY1)
		Message0("you have get the porosity %g in the cell. \n", porous_in_cell);*/
		return porous_in_cell;
	}
}

/**formula for calculate the dust cake height */
real DUST_LOAD_HEIGHT(real dust_porosity, real dust_f_load,
	real particle_den)
{

	if (dust_porosity <= 0.98 && dust_porosity >= 0.02)
	{
		return (dust_f_load / (1 - dust_porosity)) / particle_den;
	}
	else
	{

		return 0;
	}
}

int MEASSURE_TO_LAST(real* last_cord, real* meassure_cord, real distance_2, real prx_factor)
{
	real meassure_length;
	meassure_length = LENGTH_TWO_POINTS(last_cord, meassure_cord, ND_ND);
	return (COMPARE_TWO(meassure_length, distance_2, prx_factor));
}

real LENGTH_TWO_POINTS(real* one, real* two, int dim)
{
	if (dim == 3)
	{
		return(sqrt((one[2] - two[2]) * (one[2] - two[2]) + (one[1] - two[1]) * (one[1] - two[1]) + (one[0] - two[0]) * (one[0] - two[0])));
	}

	if (dim == 2)
	{
		return(sqrt((one[1] - two[1]) * (one[1] - two[1]) + (one[0] - two[0]) * (one[0] - two[0])));
	}
	else
	{
		return  (one[0] - two[0]);
	}
}

/*count the nearest cell number on the deposition face*/
int NEAR_LAYER_NUM(cell_t c, Thread* t, int deposit_face_id)
{
	int num1 = 0, n = 0;
	int tt = 0;

	begin_c_loop_int(c, t)
	{
		c_face_loop(c, t, n)
		{
			if PRINCIPAL_FACE_P(C_FACE(c, t, n), C_FACE_THREAD(c, t, n))
			{
				if (THREAD_ID(C_FACE_THREAD(c, t, n)) == deposit_face_id)
				{
					tt++;
				}
			}

		}
	}
	end_c_loop_int(c, t)


		return tt;
}

/*count the face dust load in a cell's deposition face, mainly in x y direction
* , only for 3D*/
real FACE_DUST_LOAD_C(real total_mass_c, real deposition_face)
{
	return (total_mass_c / deposition_face);
}

