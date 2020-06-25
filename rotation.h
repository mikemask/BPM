double* rotation(double *q, double *vec)
{
    double right[4];
    double q_conj[4];
    static double rot_vec[4];

    right[0] = -vec[0]*q[1] - vec[1]*q[2] - vec[2]*q[3];
    right[1] = vec[0]*q[0] + vec[1]*q[3] - vec[2]*q[2];
    right[2] = vec[1]*q[0] - vec[0]*q[3] + vec[2]*q[1];
    right[3] = vec[2]*q[0] + vec[0]*q[2] - vec[1]*q[1];

    q_conj[0] = q[0];
    q_conj[1] = -q[1];
    q_conj[2] = -q[2];
    q_conj[3] = -q[3];

    rot_vec[0] = q_conj[0]*right[0] - q_conj[1]*right[1] - q_conj[2]*right[2] - q_conj[3]*right[3];
    rot_vec[1] = q_conj[0]*right[1] + q_conj[1]*right[0] + q_conj[2]*right[3] - q_conj[3]*right[2];
    rot_vec[2] = q_conj[0]*right[2] + q_conj[2]*right[0] - q_conj[1]*right[3] + q_conj[3]*right[1];
    rot_vec[3] = q_conj[0]*right[3] + q_conj[3]*right[0] + q_conj[1]*right[2] - q_conj[2]*right[1];

    return rot_vec;
}
