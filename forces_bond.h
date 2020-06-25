void f_bonds(Bond *bond, Part parti, Part partj, double dt, double gamma)
{
    double delta_fn[3], delta_fs[3], delta_tn[3], delta_ts[3];
    double *fn, *fs, *tn, *ts;
    double *inertia;
    double *vel_i, *vel_j;
    double rel_vel[3], rel_rot[3];
    double *om_i, *om_j;
    double *pos_i, *pos_j;
    double norm[3], norm_vers[3], norm_mod, tan_vers[3], tan_mod;
    double rel_vel_n[3], rel_vel_s[3];
    double rel_rot_n[3], rel_rot_s[3];
    double R[3][3];
    double scal_vel, scal_rot, scal_forcen, scal_forces, scal_torquen, scal_torques;
    double kn, ks, radius, area;
    double theta;

//Getting bonds values

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    fn = bond -> getForce_n();
    fs = bond -> getForce_st();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

//Getting particles values and relative motions computation

    vel_i = parti.getVel();
    vel_j = partj.getVel();

    om_i = parti.getRot();
    om_j = partj.getRot();

    pos_i = parti.getPos();
    pos_j = partj.getPos();

    rel_vel[0] = vel_j[0]-vel_i[0];
    rel_vel[1] = vel_j[1]-vel_i[1];
    rel_vel[2] = vel_j[2]-vel_i[2];

    rel_rot[0] = om_j[0]-om_i[0];
    rel_rot[1] = om_j[1]-om_i[1];
    rel_rot[2] = om_j[2]-om_i[2];

//Normal versor computation

    norm_mod = distance(pos_i,pos_j);

    norm[0] = pos_j[0]-pos_i[0];
    norm[1] = pos_j[1]-pos_i[1];
    norm[2] = pos_j[2]-pos_i[2];

    norm_vers[0] = norm[0]/norm_mod;
    norm_vers[1] = norm[1]/norm_mod;
    norm_vers[2] = norm[2]/norm_mod;

//Vectors projections into normal and shear direction w.r.t. contact plane

    scal_vel = rel_vel[0]*norm_vers[0]+rel_vel[1]*norm_vers[1]+rel_vel[2]*norm_vers[2];
    scal_rot = rel_rot[0]*norm_vers[0]+rel_rot[1]*norm_vers[1]+rel_rot[2]*norm_vers[2];

    rel_vel_n[0] = scal_vel*norm_vers[0];
    rel_vel_n[1] = scal_vel*norm_vers[1];
    rel_vel_n[2] = scal_vel*norm_vers[2];

    rel_rot_n[0] = scal_rot*norm_vers[0];
    rel_rot_n[1] = scal_rot*norm_vers[1];
    rel_rot_n[2] = scal_rot*norm_vers[2];

//    rel_vel_s[0] = rel_vel[0] - rel_vel_n[0];
//    rel_vel_s[1] = rel_vel[1] - rel_vel_n[1];
//    rel_vel_s[2] = rel_vel[2] - rel_vel_n[2];
//
//    rel_rot_s[0] = rel_rot[0] - rel_rot_n[0];
//    rel_rot_s[1] = rel_rot[1] - rel_rot_n[1];
//    rel_rot_s[2] = rel_rot[2] - rel_rot_n[2];

//Computation of rotation angle and rotational matrix

    theta = sqrt((rel_rot_n[0]*rel_rot_n[0])+(rel_rot_n[1]*rel_rot_n[1])+(rel_rot_n[2]*rel_rot_n[2]))*dt;

    R[0][0] = cos(theta) + norm_vers[0]*norm_vers[0]*(1-cos(theta));
    R[0][1] = norm_vers[0]*norm_vers[1]*(1-cos(theta))-norm_vers[2]*sin(theta);
    R[0][2] = norm_vers[0]*norm_vers[2]*(1-cos(theta))+norm_vers[1]*sin(theta);

    R[1][0] = norm_vers[0]*norm_vers[1]*(1-cos(theta))+norm_vers[2]*sin(theta);
    R[1][1] = cos(theta) + norm_vers[1]*norm_vers[1]*(1-cos(theta));
    R[1][2] = norm_vers[1]*norm_vers[2]*(1-cos(theta))-norm_vers[0]*sin(theta);

    R[2][0] = norm_vers[0]*norm_vers[2]*(1-cos(theta))-norm_vers[1]*sin(theta);
    R[2][1] = norm_vers[1]*norm_vers[2]*(1-cos(theta))+norm_vers[0]*sin(theta);
    R[2][2] = cos(theta) + norm_vers[2]*norm_vers[2]*(1-cos(theta));

//Computation of rotated shear vectors

    rel_vel_s[0] = R[0][0]*(rel_vel[0]-rel_vel_n[0])+R[0][1]*(rel_vel[1]-rel_vel_n[1])+R[0][2]*(rel_vel[2]-rel_vel_n[2]);
    rel_vel_s[1] = R[1][0]*(rel_vel[0]-rel_vel_n[0])+R[1][1]*(rel_vel[1]-rel_vel_n[1])+R[1][2]*(rel_vel[2]-rel_vel_n[2]);
    rel_vel_s[2] = R[2][0]*(rel_vel[0]-rel_vel_n[0])+R[2][1]*(rel_vel[1]-rel_vel_n[1])+R[2][2]*(rel_vel[2]-rel_vel_n[2]);

    rel_rot_s[0] = R[0][0]*(rel_rot[0]-rel_rot_n[0])+R[0][1]*(rel_rot[1]-rel_rot_n[1])+R[0][2]*(rel_rot[2]-rel_rot_n[2]);
    rel_rot_s[1] = R[1][0]*(rel_rot[0]-rel_rot_n[0])+R[1][1]*(rel_rot[1]-rel_rot_n[1])+R[1][2]*(rel_rot[2]-rel_rot_n[2]);
    rel_rot_s[2] = R[2][0]*(rel_rot[0]-rel_rot_n[0])+R[2][1]*(rel_rot[1]-rel_rot_n[1])+R[2][2]*(rel_rot[2]-rel_rot_n[2]);

//Force increment computation and total force update

    delta_fn[0] = kn*area*dt*rel_vel_n[0];
    delta_fn[1] = kn*area*dt*rel_vel_n[1];
    delta_fn[2] = kn*area*dt*rel_vel_n[2];

    delta_fs[0] = -ks*area*dt*rel_vel_s[0];
    delta_fs[1] = -ks*area*dt*rel_vel_s[1];
    delta_fs[2] = -ks*area*dt*rel_vel_s[2];

    fn[0] = fn[0]*gamma + delta_fn[0];
    fn[1] = fn[1]*gamma + delta_fn[1];
    fn[2] = fn[2]*gamma + delta_fn[2];

    fs[0] = fs[0]*gamma + delta_fs[0];
    fs[1] = fs[1]*gamma + delta_fs[1];
    fs[2] = fs[2]*gamma + delta_fs[2];


    bond -> setForce_n(fn);
    bond -> setForce_st(fs);
}
