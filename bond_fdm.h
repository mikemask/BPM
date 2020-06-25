void bond_fdm(Bond *bond, Part parti, Part partj, double dt)
{
    double *pos_i, *pos_j;
    double *pos0_i, *pos0_j;
    double *q_i, *q_j;
    double kn, ks;

    double f_n[3], f_st[3], f_sr[3], t_st[3], t_sr[3], t_br[3], t_tr[3];

    double r_0[3]={}, r_c[3]={}, r_f[3]={}, gamma=0., s[3]={}, t[3]={};

//Getting bonds values

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

//Getting particles values and relative motions computation

    pos_i = parti.getPos();
    pos_j = partj.getPos();

    pos0_i = parti.getPosIn();
    pos0_j = partj.getPosIn();

    q_i = parti.getQuat();
    q_j = partj.getQuat();

    r_f[0] = pos_i[0] - pos_j[0];
    r_f[1] = pos_i[1] - pos_j[1];
    r_f[2] = pos_i[2] - pos_j[2];

    r_0[0] = pos0_i[0] - pos0_j[0];
    r_0[1] = pos0_i[1] - pos0_j[1];
    r_0[2] = pos0_i[2] - pos0_j[2];

    //double *q_rc;

    //q_rc = rotation(q_j, r_f);

    r_c[0] = r_f[0];
    r_c[1] = r_f[1];
    r_c[2] = r_f[2];

    double abs_rc = abs(r_c);
    double abs_r0 = abs(r_0);

    f_n[0] = kn*(abs_rc - abs_r0)/abs_rc*r_c[0];
    f_n[1] = kn*(abs_rc - abs_r0)/abs_rc*r_c[1];
    f_n[2] = kn*(abs_rc - abs_r0)/abs_rc*r_c[2];

    bond -> setForce_n(f_n);

    gamma = acos((r_0[0]*r_c[0]+r_0[1]*r_c[1]+r_0[2]*r_c[2])/(abs_rc*abs_r0));

    double S[3]={};

    S[0] = r_c[1]*(r_c[0]*r_0[1]-r_0[0]*r_c[1]) - r_c[2]*(r_0[0]*r_c[2]-r_c[0]*r_0[2]);
    S[1] = r_c[2]*(r_c[1]*r_0[2]-r_0[1]*r_c[2]) - r_c[0]*(r_0[1]*r_c[0]-r_c[1]*r_0[0]);
    S[2] = r_c[0]*(r_c[2]*r_0[0]-r_0[2]*r_c[0]) - r_c[1]*(r_0[2]*r_c[1]-r_c[2]*r_0[1]);

    double abs_S = abs(S);

    if (abs_S < 1.e-8)
    {
        f_st[0] = 0.;
        f_st[1] = 0.;
        f_st[2] = 0.;
    }

    else
    {
        s[0] = S[0]/abs_S;
        s[1] = S[1]/abs_S;
        s[2] = S[2]/abs_S;

        f_st[0] = ks*abs_r0*gamma*s[0];
        f_st[1] = ks*abs_r0*gamma*s[1];
        f_st[2] = ks*abs_r0*gamma*s[2];
    }

    bond -> setForce_st(f_st);

}

